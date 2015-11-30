#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
    make a new BEAST  config file by inserting FASTA sequences from the RV217 trail.
    The sequence names int the RV217 data look like this:
	>RV217_PDB|1M|01WG|NFLG|2011/11/10

    The 5th component is the sample data and should eb converted to an appropriate tipdate for BEAST.

    usage:

         mkbeast.py -p foo template.xml sequences.fasta  >beast_in.xml

    This will take the fasta sequences and insert them into template.xml
    to produce an XML file that is suitable to pass to BEAST.  

    Once the BEAST config file is generated, you would run,

         beast beast_in.xml

    This will produce various output files, among which is foo.trees.
    That file gets fed to the 'annotatetrees' program and the output of that gets
    visualized with 'figtree'.
'''

from __future__ import print_function

import sys
import os.path
import argparse
from Bio import SeqIO

def deduplicate(datafile):
    '''
    Read sequences from a FASTA file (datafile) and create a nested data structure thet organizaes the sequences by patient and sample date.
    '''
    seqs = dict()
    counts = dict()
    with open(datafile, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta") :
            h = hash(str(record.seq))
            if h not in seqs:
                seqs[h] = record
                counts[h] = 1
            else:
                counts[h] += 1
    for k,r in seqs.items():
        r.id = "{}_{}".format(r.id, counts[k])
        r.description = ""
        yield(r)

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

    parser.add_argument('datafiles', nargs='+', help='FASTA input', type=existing_file)

    return parser


def main(args=sys.argv[1:]):

    parser = build_parser()
    a = parser.parse_args()
    
    for datafile in a.datafiles:
        SeqIO.write(deduplicate(datafile), sys.stdout, "fasta")

    
if __name__ == "__main__":
   main(sys.argv[1:])
   

