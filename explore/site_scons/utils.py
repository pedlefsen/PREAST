

# read in the csv
# pick out specific columns
# skip burnin period
# calculate mean.

import sys
import argparse
import csv
import numpy as np

from SCons.Action import ActionFactory
import SCons.Util
import os
import time

def parse_log(fp, burnin=0.9):
    """
    generator function for parsing ancestral founder sequences from BEAST trait file.

    The file expected to be a set of ancestral founder sequences from Beast.
    Each non-comment line consists of an iteration count followed by some number or quoted founder sequences.
    """

    # first count how many sequences are in the file so we can figure out how many to skip for burn-in.
    count = 0
    for line in csv.reader(filter(lambda row: row[0]!='#', fp), delimiter='\t'):
        count += 1

    skip = round(count * burnin)
    
    # Reposition pointer at the beginning once again
    fp.seek(0, 0);
    for line in csv.DictReader(filter(lambda row: row[0]!='#', fp), delimiter='\t'):
        if skip > 0:
            skip -= 1
            continue
        yield line


def logfile_iter(fh, column, burnin):
    """
    iterate over BEAST log file

    This routine knows how to retrieve and iterate over logfile
    entries.

    :param path: string name of directory where sequence files can be found.
    :returns: iterator over values
    :rtype:
    """

    # NB parse_log() will skip a burn-in period
    for row in parse_log(fh, burnin=burnin):
        yield row[column]

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--burnin', type=int, help='skip this many entries as a burn-in period (default: no burn-in period)',
            action='store', default=0, dest='burnin')
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin)
    return parser

def extract_posterior(path, column='clock.rate', burnin=0.4):
    with open(str(path), 'rb') as fp:
        generator = logfile_iter(fp, column, burnin)
        values = np.array(list(generator)).astype(np.float)
    return np.mean(values)



# Running commands on the cluster sometimes has the unfortunate side-effect of
# letting distributed filesystems get out of sync.  A file that is written on
# the cluster may not be visible on local machines for several seconds.  This
# doesn't happen all the time, and it can be difficult to demonstrate, but it
# often occurs when local and cluster commands manipulate the same file in
# quick succession.  The Wait() action is meant to wait for a file to become
# locally visible after running a command on the cluster (via srun or salloc).
#
# Wait() usage is typically,
# 	target='output.file'
# 	env.Command(target, 'source.file',
#   	        [ "srun some-command <${TARGET}",
#    			   Wait(target)
#				])
#
# This will cause the execution to pause after running 'some-command' until the target shows up on the local machine.
# The target will be polled on a 2-second interval, and the command will fail if the target does not show up within about 10 seconds.




def get_paths_str(dest):
    # If dest is a list, we need to manually call str() on each element
    if SCons.Util.is_List(dest):
        elem_strs = []
        for element in dest:
            elem_strs.append('"' + str(element) + '"')
        return '[' + ', '.join(elem_strs) + ']'
    else:
        return '"' + str(dest) + '"'

# https://github.com/azatoth/scons/blob/73f996e59902d03ec432cc662252aac5fb72f1f8/src/engine/SCons/Defaults.py 
def wait_func(dest):
    SCons.Node.FS.invalidate_node_memos(dest)
    if not SCons.Util.is_List(dest):
        dest = [dest]
    for entry in dest:
        count = 0
        limit = 3
        while not os.path.isfile(entry) or os.stat(entry).st_size == 0:
            print("waiting for {}...".format(entry))
            time.sleep(2)
            count = count + 1
            if count >limit:
                return 1
    return 0

Wait = ActionFactory(wait_func, lambda dir: 'Wait(%s)' % get_paths_str(dir))

            


            
