#!/usr/bin/env bash
set -e # exit on first error
if [ ! -e venv/bin/python2.7 ]; then
    virtualenv -q venv
    . venv/bin/activate
    pip -q install -r requirements.txt
fi
echo 'module use modules;'
echo 'module load BEAST/1.8.2 beagle/2.1;'
echo 'source venv/bin/activate;'

