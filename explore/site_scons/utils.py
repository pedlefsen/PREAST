#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from SCons.Action import ActionFactory
import SCons.Util
import os
import time


def extract_posterior(path, column='clock.rate', burnin=0.4):
    df = pd.read_csv(path, sep='\t', comment='#')
    return df[column][int(len(df.index)*burnin):].mean()


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

            


            
