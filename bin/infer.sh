#!/usr/bin/env bash
set -e # exit on error

# see http://unix.stackexchange.com/a/84980/100709
# for a mktmpdir that works on both linux and OS X
outdir=$(mktemp -d 2>/dev/null)
function cleanup() {
    if [ -d "${outdir}" -a "${outdir}"=="$(dirname $(mktemp -d -u))*" ] && [ ${keep} != true ]
    then
	echocmd rm -rf "${outdir}"
    fi
}

trap "cleanup" EXIT

export PATH=~matsengrp/local/bin/:$PATH
export rate=1.62e-2
export backoff=0

function usage() {
cat << EOF
Usage: ${0##*/} [-v] [-d] [-p <prefix>] <sequences.fasta>
Infer founders and timing for HIV samples.

    -h          display this help and exit
    -v		verbose mode.  Echo description of each command being run.
    -k		keep temporary files. Useful in conjunction with -v to debug the process.
    -p <prefix>:	supply a prefix to use for output files (default: "infer_").
EOF

}

function echotty() {
    # echo to stderr if bound to a tty, otherwise stay silent.
    if [  ${verbose} = true ] && [ -t 2 ] ; then
	echo "$@" 1>&2
    fi
}

function echocmd() {
    if [ ${verbose} = true ]; then
	echo "$@" 1>&2
    fi
    eval $@
}


# get full directory name of the script
# http://stackoverflow.com/a/246128/1135316
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )"  >/dev/null && pwd )"

export PATH=${DIR}:${PATH}
export TEMPLATES=${DIR}/../templates
files=()

# Initialize our own variables:
output_file=""
verbose=false
keep=false
prefix="infer_"

OPTIND=1 # Reset is necessary if getopts was used previously in the script.  It is a good idea to make this local in a function.
while getopts "hvkp:" opt; do
    case "$opt" in
        h)
            usage
            exit 0
            ;;
        v)  verbose=true
            ;;
        k)  keep=true
            ;;
        p)  prefix="$OPTARG"
            ;;
        '?')
            usage >&2
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))" # Shift off the options and optional --.

files=$@

if [[ ${#files[@]} == 0 ]]; then
    usage >&2
    exit 1
fi

mkdir -p $(dirname ${prefix})

for sample in "${files[@]}"
do
    label=$(basename $sample)
    label=${label%.*}
    echocmd "dedup.py ${sample} >${outdir}/${label}.fa"

    echotty "Creating alignment with PRANK..." 
    echocmd "prank -d=${outdir}/${label}.fa -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"

    # multiple founders
    # split just below the root and sequences separately.
    echocmd "treesplit.py -o  ${outdir} ${outdir}/prank.best.dnd ${outdir}/${label}.fa"
    
    if [[ $(egrep -c '^>' ${outdir}/left.fasta) > 1 ]]
    then
	echocmd "prank -d=${outdir}/left.fasta -o=${outdir}/left -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"
	founder1=$(prankroot.py ${outdir}/left.best.anc.dnd ${outdir}/left.best.anc.fas)
    else
	founder1=$(egrep -v '^>' ${outdir}/left.fasta)
    fi
    
    if [[ $(egrep -c '^>' ${outdir}/right.fasta) > 1 ]]
    then
	echocmd "prank -d=${outdir}/right.fasta -o=${outdir}/right -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"
	founder2=$(prankroot.py ${outdir}/right.best.anc.dnd ${outdir}/right.best.anc.fas)
    else
	founder2=$(egrep -v '^>' ${outdir}/right.fasta)
    fi
    
    # if [[ "${founder1}" == "${founder2}" ]]
    # then
    #     echo "Founders are identical"
    # else
    #     echo "Founders are different"
    # fi
    cat >${prefix}multiple.fa <<EOF
>${label}_founder1
${founder1}
>${label}_founder2
${founder2}
EOF
    
    # single founder - process all sequences together
    # extract sequences at root of guide tree and output in FASTA format.
    cat >${prefix}single.fa <<EOF
>${label}_founder
`prankroot.py ${outdir}/prank.best.anc.dnd ${outdir}/prank.best.anc.fas`
EOF
    
	
    # inserts sample sequences into a beast configurarion file.
    # parses date fromt he sequence ids and converts to tip date for BEAST
    echotty "Make BEAST config file..."
    echocmd "mkbeast_rv217.py --template ${TEMPLATES}/beast_strict.template rate=${rate} backoff=${backoff} ${outdir}/prank.best.fas > ${outdir}/beast_in.xml"

    # sigh...if using the -working option of beast, the configuration file must be specified with an absoute path.
    echotty "Running BEAST..."
    echocmd "beast -working -overwrite -beagle ${outdir}/beast_in.xml >${outdir}/beastcmd.log  2>&1"

    echotty 'Extracting estimated time of infection'
    echocmd 'posterior_toi.py ${outdir}/beastout.log ${outdir}/${label}.fa' >${prefix}toi.csv

    # echotty "Exracting guide tree..." 
    # echocmd "treeannotator ${outdir}/beastout.trees > ${outdir}/mcc.tree 2>${outdir}/treeannotator.log"

    # tr "/" "-" < ${outdir}/mcc.tree | nexus2newick.py | tr "-" "/" | tr -d "'" > ${outdir}/guidetree.tree

    # echotty "Inferring ancestral states with PRANK..." 
    # echocmd "prank -d=${outdir}/${label}.fa -t=${outdir}/guidetree.tree -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -codon >${outdir}/prankcmd.log 2>&1"

done

