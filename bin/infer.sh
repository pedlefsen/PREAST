#!/usr/bin/env bash
set -e # exit on error

# see http://unix.stackexchange.com/a/84980/100709
# for a mktmpdir that works on both linux and OS X
outdir=$(mktemp -d 2>/dev/null)
function cleanup() {
    if [ -d "${outdir}" -a "${outdir}"=="$(dirname $(mktemp -d -u))*" ] && [ ${debug} != true ]
    then
	echocmd rm -rf "${outdir}"
    fi
}

trap "cleanup" EXIT

export PATH=~matsengrp/local/bin/:$PATH



function usage() {
    echo "infer.sh <sequences.fasta>";
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

debug=false
verbose=false

while [[ $# > 0 ]]
do
key="$1"
case $key in
    -v|--verbose)
	verbose=true
    ;;
    -d|--debug)
	debug=true
    ;;
    -*)
	echo 'Unknown argument - "$key"'
	usage
	exit 1
	;;
    *)
	files+=($key)
    ;;
esac
shift # past argument
done


if [[ ${#files[@]} == 0 ]]; then
    usage
    exit 1
fi

for sample in "${files[@]}"
do
    label=$(basename $sample)
    label=${label%.*}
    echocmd "dedup.py ${sample} >${outdir}/${label}.fa"

    echotty "Creating alignment with PRANK..." 
    echocmd "prank -d=${outdir}/${label}.fa -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"

    echotty "Testing for multiple founders..."
    echocmd "keele ${outdir}/${label}.fa >${outdir}/keele.out"
    moi=$(cat ${outdir}/keele.out)
    if [[ "${moi}" == *"multiple infection"* ]]
    then
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
	echo ">" ${label}_founder_1
	echo ${founder1}
	echo ">" ${label}_founder_2
	echo ${founder2}

    else
	# single founder - process all sequences together
	# extract sequences at root of guide tree and output in FASTA format.
	echo ">" ${label}_founder
	echocmd "prankroot.py ${outdir}/prank.best.anc.dnd ${outdir}/prank.best.anc.fas"
    fi
    
	
    # inserts sample sequences into a beast configurarion file.
    # parses date fromt he sequence ids and converts to tip date for BEAST
    echotty "Make BEAST config file..."
    echocmd "mkbeast_rv217.py --template ${TEMPLATES}/beast_strict.template ${outdir}/prank.best.fas > ${outdir}/beast_in.xml"

    # sigh...if using the -working option of beast, the configuration file must be specified with an absoute path.
    echotty "Running BEAST..."
    echocmd "beast -working -overwrite -beagle ${outdir}/beast_in.xml >${outdir}/beastcmd.log  2>&1"

    echotty 'Extracting estimated time of infection'
    echocmd 'posterior_toi.py ${outdir}/beastout.log ${outdir}/${label}.fa'

    # echotty "Exracting guide tree..." 
    # echocmd "treeannotator ${outdir}/beastout.trees > ${outdir}/mcc.tree 2>${outdir}/treeannotator.log"

    # tr "/" "-" < ${outdir}/mcc.tree | nexus2newick.py | tr "-" "/" | tr -d "'" > ${outdir}/guidetree.tree

    # echotty "Inferring ancestral states with PRANK..." 
    # echocmd "prank -d=${outdir}/${label}.fa -t=${outdir}/guidetree.tree -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -codon >${outdir}/prankcmd.log 2>&1"

done

