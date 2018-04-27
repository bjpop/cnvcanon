#!/usr/bin/env python

from distutils.core import setup

LONG_DESCRIPTION = \
'''XXX fixme'''


setup(
    name='cnvcanon',
    version='0.1.0.0',
    author='Bernie Pope',
    author_email='bjpope@unimelb.edu.au',
    packages=['cnvcanon'],
    package_dir={'cnvcanon': 'cnvcanon'},
    entry_points={
        'console_scripts': ['cnvcanon = cnvcanon.cnvcanon:main',
            ]
    },
    url='https://github.com/bjpop/cnvcanon',
    license='LICENSE',
    description=('XXX fixme'),
    long_description=(LONG_DESCRIPTION),
    install_requires=["networkx", "quicksect"],
)
