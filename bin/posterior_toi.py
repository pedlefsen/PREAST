#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
extract an estimate of time of infection from a mcmc posterior distribution.

Note at one time this script parsed the FASTA file to produce an actual
date of infection.  Now it just produces treeheight in days, so no
fasta file is necessary.

Typical usage:

$ posterior_toi.py outout/beastout.log  >toi.csv
'''
from __future__ import print_function

import re
from Bio import SeqIO
from collections import defaultdict
from datetime import datetime, date, timedelta
import calendar
import sys
import argparse
import os.path
import pandas as pd

# Convert a string like '1M', '3M', or '6M' to a timedelta object corresponding to number of months indicated in string.
#
def str2timedelta(s):
    intervals = { 'D' : 1,
                  'W' : 7,
                  'M' : 30,
                  'Y' : 365 }

    m = re.match(r'(\d+)([MWYD])$', s)
    if m.group(2) in intervals:
        delay = int(m.group(1)) * intervals[m.group(2)]
    else:
        delay = 30 * 12
    return(timedelta(days=delay))

    

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

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--burnin', '-b', help='burning percentage, expressed as whole integer [default %default]', default=20)
    parser.add_argument('logfile', nargs=1, help='Beast log file', type=existing_file)

    return parser


def main(args=sys.argv[1:]):
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()
    
    tbl = pd.read_table(a.logfile[0], skiprows=2)
    rootHeight = tbl['treeModel.rootHeight']
    n = len(rootHeight)
    burnin = int(n*(a.burnin / 100.0))
    toi = tbl['treeModel.rootHeight'][burnin:].mean()
    std = tbl['treeModel.rootHeight'][burnin:].std()


    # calculate a 95% credible intervale (+- 1.96 stddev),
    # assuming treeheights are normally distributed.
    # doi = latest_timepoint - timedelta(days=365.0*toi)
    # doi_early = latest_timepoint - timedelta(days=365.0*(toi + 1.96*std))
    # doi_late = latest_timepoint - timedelta(days=365.0*(toi - 1.96*std))

    print('{:.0f},{:.0f},{:.0f}'.format(
        365.0*toi, 365.0*(toi + 1.96*std), 365.0*(toi - 1.96*std)))
        #print(tbl)

    
if __name__ == "__main__":
   main(sys.argv[1:])
   
