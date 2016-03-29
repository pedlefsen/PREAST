#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
extract-params - extract mean posterior values after burnin from BEAST log files.
				 Output those values as JSON.

The purpose of this script is to extract values from the posterior of
a BEAST run that will be used to configure a subsequent BEAST run.
"""
from __future__ import print_function

import sys
from argparse import ArgumentParser, RawTextHelpFormatter
import pandas as pd

def extract_posterior(path, column='clock.rate', burnin=0.4):
    df = pd.read_csv(path, sep='\t', comment='#')
    return df[column][int(len(df.index)*burnin):].mean()

def build_parser():
    """
    Build the command-line argument parser.
    """

    def existing_file(fname):
        """
        Argparse type for an existing file
        """
        if not os.path.isfile(fname):
            raise ValueError("Invalid file: " + str(fname))
        return fname

    parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
    parser.add_argument('files', nargs=1, help='BEAST statistics log file, e.g. beastout.log')
    parser.add_argument('cols', nargs='*', default=['clock.rate', 'exponential.growthRate'],
                            help="Columns to extract as space-separated list,\n"
                                  "e.g. 'clock.rate' 'exponential.growthRate'")

    return parser

def main(args=sys.argv[1:]):
    parser = build_parser()
    a = parser.parse_args()

    burnin=0.4
    df = pd.read_csv(a.files[0], sep='\t', comment='#')
    params = df[a.cols][int(len(df.index)*burnin):].mean()
    print(params.to_json())
    
if __name__ == "__main__":
   main(sys.argv[1:])
