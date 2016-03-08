#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
make a new BEAST  config file by inserting FASTA sequences from the RV217 trial.

The sequence names in the RV217 data look like this:
>RV217_PDB|1M|01WG|NFLG|2011/11/10

The 5th component is the sample data and should be converted to an appropriate tipdate for BEAST.

usage:

mkbeast.py -p foo template.xml sequences.fasta  >beast_in.xml

This will take the fasta sequences and insert them into template.xml
to produce an XML file that is suitable to pass to BEAST.  

Once the BEAST config file is generated, you would run,

beast beast_in.xml

This will produce various output files, among which is foo.trees.
That file gets fed to the 'annotatetrees' program and the output of that gets
visualized with 'figtree'.

Templates are Beast XML config files with embedded Jinja2 fragments
that are expanded and substituted before rendering.
'''

from __future__ import print_function

import jinja2

import re
from Bio import SeqIO
from collections import defaultdict
import getpass
from datetime import datetime, date, timedelta
import calendar
import sys
import argparse
import os.path
import operator

# patients dict stores instance of class Patient keyed by patient
# id string.  we create a new one when we run into a patient id
# string that we haven't seem before.  patient id string is
# derived from the first part of the name on an alignment,
# e.g. 'patient1_100_1' yields patient id 'patient1'
#
def dayofyear(value):
    t = value.timetuple()
    return "{}".format(t.tm_year+(t.tm_yday-1)/(366.0 if calendar.isleap(t.tm_year) else 365.0))

def dateformat(value, format='%d/%m/%Y'):
    return value.strftime(format)

def deltayears(value):
    return value.total_seconds()/(60*60*24*365)

def render(params, template, fp):
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(searchpath="/"), undefined=jinja2.StrictUndefined)
    env.filters['dayofyear'] = dayofyear
    env.filters['deltayears'] = deltayears
    env.filters['dateformat'] = dateformat
    
    # Alias str.format to strformat in template
    env.filters['strformat'] = str.format
    template = env.get_template(os.path.abspath(template))

    # define additional values available in the template
    vars = dict(date=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                  user=getpass.getuser(),
                  command=" ".join(sys.argv),
                  workdir=os.getcwd(),
                  now=datetime.now(),
                  timedelta=timedelta)
    # merge in the values passed in to this routine.
    vars.update(params)
    
    template.stream(**vars).dump(fp)

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

sample_dd = lambda: dict()
patient_dd = lambda: defaultdict(sample_dd)
patients_dd = lambda: defaultdict(patient_dd)

byseq_dd = lambda: defaultdict(set)
bydate_dd = lambda: defaultdict(byseq_dd)


def processFasta(file):
    bydate = bydate_dd()
    with open(file, "rU") as fh:
        # for each fasta sequence in the data file, create a taxon node and a sequence node.
        for r in SeqIO.parse(fh, "fasta"):
            fields = r.id.split('|')
            patientId = fields[0]
            timePoint = fields[1] if len(fields) > 0 else "0"
            sampleDate = fields[4] if len(fields) > 3 else "0"
            # if the sequence name ends with "_<digits>", assume this is a multiplicity indicator
            # and strip it off before converting the date.
            sampleDate = re.sub(r"_\d+$", '', sampleDate)
            sampleDate = datetime.strptime(sampleDate, '%Y/%m/%d')

            byseq = bydate[sampleDate]
            byseq[str(r.seq)].add(str(r.id))

    if not bydate:
        raise Exception('Empty fasta file - {}'.format(file))

    return bydate


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

    parser.add_argument('--template', '-t', help='templated BEAST XML config',
            required=True, dest='template')
    parser.add_argument('--prefix', help='Specify a prefix for all output log filename',
            default="", dest='prefix')
    parser.add_argument('--mono', '-m', help='Duplicate sequences from the sametimepoint are marked as monophyletic', action='store_true',
            default=False, dest='mono')
    parser.add_argument('--nodata', action='store_true', help='set all sequences to "N".  Useful for exploring prior distributions without data.')
    class StoreNameValuePair(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            keyvalues = [ v for v in values if ('=' in v and not os.path.isfile(v)) ]
            values = [ v for v in values if not v in keyvalues]
            setattr(namespace, 'keyvalues', dict([kv.split('=') for kv in keyvalues ]))
            setattr(namespace, self.dest, values)

    
    parser.add_argument('file', nargs='+', help='FASTA input', action=StoreNameValuePair)

    return parser


def main(args=sys.argv[1:]):
    '''
    Parse a generic template and insert sequences from a FASTA file into the middle,
    separated by the appropriate XML tags.
    '''

    parser = build_parser()
    a = parser.parse_args()

    if 'toi' in a.keyvalues:
        toi = a.keyvalues['toi']
        toi = toi.split(',')
        if len(toi) != 2:
            print('Bad toi value - expected two integers values separated by comma', file=sys.stderr);
            sys.exit(1)
        toi = [int(t) for t in toi]
        if sum([t < 0 for t in toi]) > 0:
            print('Bad toi value - toi priors must be strictly positive values', file=sys.stderr);
            sys.exit(1)
        if toi[0] >= toi[1] :
            print('Bad toi value - first toi value must be less than second', file=sys.stderr);
            sys.exit(1)
        a.keyvalues['toi'] = [timedelta(days=int(days)) for days in toi]


    bydate = processFasta(a.file[0])

    earliest_timepoint = min(bydate.keys())
    latest_timepoint = max(bydate.keys())
    
    params = dict(sequences=bydate,
                      earliest_timepoint = earliest_timepoint,
                      latest_timepoint = latest_timepoint,
                      monophyly=a.mono)
    params.update(a.keyvalues)

    render(params, a.template, sys.stdout)



    
if __name__ == "__main__":
   main(sys.argv[1:])
   


