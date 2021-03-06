'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv
from quicksect import IntervalTree
import networkx as nx
from pathlib import Path
from copy import copy
import re
import plotly as py
import plotly.graph_objs as go
from plotly.graph_objs import Scatter, Layout


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_CSV_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_OVERLAP = 0.75 
PROGRAM_NAME = "cnvcanon"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Merge distilled SVs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--overlap',
                        metavar='PROPORTION',
                        default=DEFAULT_OVERLAP,
                        type=float,
                        help='proportion overlap for CNV equality (default {})'.format(DEFAULT_OVERLAP))
    parser.add_argument('--heatmap',
                        metavar='HTML_FILE',
                        type=str,
                        help='Output a heatmap file with this name')
    parser.add_argument('--histogram',
                        metavar='HTML_FILE',
                        type=str,
                        help='Output a histogram of confidence scores with this name')
    parser.add_argument('csv_file',
                        metavar='CSV_FILE',
                        type=str,
                        help='Input CSV files')
    return parser.parse_args()


class CNVIntervals(object):
    def __init__(self):
        self.chroms = {}

    def add(self, chrom, start, end, val):
        if chrom not in self.chroms:
            self.chroms[chrom] = IntervalTree()
        tree = self.chroms[chrom]
        tree.add(start, end, val)

    def lookup(self, chrom, start, end):
        if chrom in self.chroms:
            return self.chroms[chrom].search(start, end)
        else:
            return [] 


# mapping from unique integer (count) to variant record
class Variants(object):
    def __init__(self):
        self.variants = {}
        self.count = 0

    def add(self, variant):
        self.variants[self.count] = variant
        self.count += 1


def read_csv_file(options):
    variants = Variants()
    sample_ids = set()
    with open(options.csv_file) as file:
        logging.info("Processing TSV file from %s...", options.csv_file)
        reader = csv.DictReader(file)
        for row in reader:
            variants.add(row)
            sample_ids.add(row['sampleID'])
        logging.info("Processing TSV file from %s: done", options.csv_file)
    return sample_ids, variants


def cnv_intervals(variants):
    logging.info("Computing %i CNV intervals", len(variants))
    intervals = CNVIntervals()
    for variant_id, variant_info in variants.items():
        chrom = variant_info['chr']
        start = int(variant_info['start'])
        end = int(variant_info['end'])
        intervals.add(chrom, start, end, variant_id)
    logging.info("Computing %i CNV intervals, done", len(variants))
    return intervals


def is_overlap(start1, end1, start2, end2, min_overlap):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        overlap_size = float((overlap_end - overlap_start) + 1)
        cnv1_size = (end1 - start1) + 1
        cnv2_size = (end2 - start2) + 1
        cnv1_overlap = overlap_size / cnv1_size
        cnv2_overlap = overlap_size / cnv2_size
        return cnv1_overlap >= min_overlap and cnv2_overlap >= min_overlap
    return False

def get_intersections(overlap, variants, intervals):
    logging.info("Computing %i CNV intersections...", len(variants))
    overlaps = nx.Graph() 
    for variant_id, variant_info in variants.items():
        # make sure all variants are recorded in the graph
        overlaps.add_node(variant_id)
        chrom = variant_info['chr']
        start = int(variant_info['start'])
        end = int(variant_info['end'])
        this_cn = variant_info['cn'] 
        this_conf = variant_info['conf'] 
        intersections = set(intervals.lookup(chrom, start, end))
        for other_variant in intersections:
            other_variant_id = other_variant.data
            other_variant_info = variants[other_variant_id]
            # don't add self edges
            if variant_id != other_variant_id and \
               is_overlap(start, end, other_variant.start, other_variant.end, overlap) and \
               this_cn == other_variant_info['cn']:
                overlaps.add_edge(variant_id, other_variant_id)
    logging.info("Computing %i variant intersections: done", len(variants))
    return overlaps


def list_median(items):
    mid_pos = len(items) // 2
    return sorted(items)[mid_pos] 

def average(items):
    return (sum(items) / len(items))


def merge_overlaps(sample_ids, variants, overlaps):
    logging.info("Merging overlapping variants...")
    sorted_sample_ids = sorted(sample_ids)
    writer = csv.writer(sys.stdout, delimiter="\t")
    header = ["chr", "start", "end", "cn", "state", "num_samples"] + sorted_sample_ids
    writer.writerow(header)
    results = []
    samples_cnvs = {}
    variant_ids = set()
    for component in nx.connected_components(overlaps):
        if len(component) > 0:
            variant_infos = [variants[id] for id in component]
            samples_confs = {var['sampleID']: var['conf'] for var in variant_infos}
            first_info = variant_infos[0]
            chrom  = first_info['chr']
            cn = first_info['cn']
            state = first_info['state']
            start = min([int(info['start']) for info in variant_infos])
            end = max([int(info['end']) for info in variant_infos])
            results.append([chrom, start, end, cn, state, len(component)] + [samples_confs])
    results.sort(key = lambda x: (int(x[0][3:]), int(x[1]), int(x[2]), int(x[3])))
    samples_confs = [] 
    for row in results:
        cnv_info = row[:6]
        samples_info = row[6]
        sample_output = [] 
        for sample_id in sorted_sample_ids: 
            sample_output.append(samples_info.get(sample_id, ''))
        writer.writerow(cnv_info + sample_output)
        samples_confs.append(sample_output)
    logging.info("Merging overlapping variants: done")
    return samples_confs

def plot_results(heatmap, histogram, sample_ids, samples_confs):
    sorted_sample_ids = sorted(sample_ids)
    heatmap_rows = []
    non_zeros = []
    for row in samples_confs:
        this_row = [float(conf) if conf else 0 for conf in row]
        heatmap_rows.append(this_row)
        for val in this_row:
            if val > 0:
                non_zeros.append(val)

    variant_nums = [count for count,_ in enumerate(heatmap_rows[0])]

    maximum_conf = max(non_zeros)
    tick = 0
    tick_vals = []
    tick_step = 50
    while tick <= maximum_conf:
        tick_vals.append(tick)
        tick += tick_step

    heatmap_data = [{
        'x': sorted_sample_ids,
        'y': variant_nums, 
        'z': heatmap_rows,
        'type': 'heatmap',
        'colorscale': [
            [0, 'rgb(250, 250, 250)'],        #0
            [1./10000, 'rgb(200, 200, 200)'], #10
            [1./1000, 'rgb(150, 150, 150)'],  #100
            [1./100, 'rgb(100, 100, 100)'],   #1000
            [1./10, 'rgb(50, 50, 50)'],       #10000
            [1., 'rgb(0, 0, 0)'],             #100000

        ],
        'colorbar': {
            'tick0': 0,
            'tickmode': 'array',
            'tickvals': tick_vals 
        }
    }]

    heatmap_layout = {'title': 'sample vs confidence heatmap',
        'xaxis': {'title': 'sample'}, 'yaxis': {'title': 'variant'}}
    heatmap_fig = {'data': heatmap_data, 'layout': heatmap_layout}
    py.offline.plot(heatmap_fig, filename=heatmap)

    histogram_data = [{
        'x' : non_zeros,
        'type': 'histogram'
    }]
    histogram_layout = {'title': 'histogram of PennCNV confidence scores',
        'xaxis': {'title': 'confidence scores'}, 'yaxis': {'title': 'count'}}
    histogram_fig = {'data': histogram_data, 'layout': histogram_layout}
    py.offline.plot(histogram_fig, filename=histogram)


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is None:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
    else:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('computation started')
    logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    sample_ids, variants = read_csv_file(options)
    intervals = cnv_intervals(variants.variants)
    overlaps = get_intersections(options.overlap, variants.variants, intervals)
    samples_confs = merge_overlaps(sample_ids, variants.variants, overlaps)
    plot_results(options.heatmap, options.histogram, sample_ids, samples_confs)
    logging.info("computation ended")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
