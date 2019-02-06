#!/usr/bin/env bash

# Initialize our own variables:
VERSION="0.16"
output_file=""
verbose=false
keep=false
deduplicate=false
prefix="infer_"
dryrun=false

# these are default parameters used if the user deos not supply a -params option
declare -a params
params=( ["clock_rate"]="0.0039439059" ["exponential_growthRate"]="3.6741292922")

# see http://unix.stackexchange.com/a/84980/100709
# for a mktmpdir that works on both linux and OS X
outdir=$(mktemp -d 2>/dev/null || mktemp -d -t 'mytmpdir')
function cleanup() {
    if [ -d "${outdir}" -a "${outdir}"=="$(dirname $(mktemp -d -u -t 'mytmpdir'))*" ] && [ ${keep} != true ]
    then
	echocmd rm -rf "${outdir}"
    fi
}

#trap "cleanup" EXIT

export PATH=~matsengrp/local/bin/:$PATH

function usage() {
cat << EOF
Usage: ${0##*/} [options] <sequences.fasta>
Infer founders and timing for HIV samples.

Options are,
    -h,--help		display this help and exit
    -v,--verbose	toggle verbose mode.  Echo description of each command being run.
    -n,--dryrun		show what would happen, wihout actually doing anything. Implies --verbose.
    -k,--keep		keep temporary files. Useful in conjunction with -v to debug the process.
    -d,--deduplicate	collapse duplicate sequence.

    -p <prefix>,
    --prefix <prefix>:	output file prefix, default='_infer',

    -j <params>,
    --json <params>:	supply a json parameter file,

    -t <toi>,
    --toi <toi>:	supply prior limits on time of infection (default: none).
			'toi' is a comma-seperates pair of integers representing the 
			lower- and upper-bounds on time of infection (in days) prior to 
			the most recent sample date as specified in the sequence names.

Output:
    Produces a 'infer_toi.csv' with an estimate of time of infection
    measured in dayssince since latest sample date (aka the height of
    the inferred geneaology).  Alos produces 'infer_single.fa' and
    'infer_multiple.fa' which are reconstructions of the founder
    viruses for putative single or multiple founders.

Examples:
    $ bin/infer.sh -v  bakeoff_analysis_rv217_40094_6M_NFLG_removeRecombinedSequences.fasta  

    $ bin/infer.sh -v  --toi 150,200 bakeoff_analysis_rv217_40094_6M_NFLG_removeRecombinedSequences.fasta  
EOF

}

function echocmd() {
    if [[ ${verbose} == true ]]; then
	echo "$@" 1>&2
    fi
    if [[ ${dryrun} != true ]]; then
	eval $@
    fi
}

function echotty() {
    # echo to stderr if bound to a tty, otherwise stay silent.
    if [[ -t 2 ]] ; then
	echocmd "echo $@" 
    fi
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

# check for updates.
# if checksum of this script and the version in the official repository
# do not matchs, print a warning.
reposum=$(curl --silent  https://raw.githubusercontent.com/matsengrp/PREAST/master/bin/infer.sh | md5deep -q)
thissum=$(cat $DIR/infer.sh | md5deep -q)
if [[ $reposum != $thissum ]]
then
    echo "** An update to infer.sh is avaiable at https://github.com/matsengrp/PREAST.git **"
fi

export PATH=${DIR}:${PATH}
export TEMPLATES=${DIR}/../templates
files=()

# using getopt -- http://stackoverflow.com/a/7948533/1135316
TEMP=`getopt -o hvkdp:t:j:n --long help,verbose,keep,deduplicate,prefix:,toi:,version,dryrun \
             -n 'infer' -- "$@"`

if [ $? != 0 ] ; then usage >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true; do
  case "$1" in
    -h | --help ) usage; exit 0 ;;
    -v | --verbose ) verbose=true; shift ;;
    -k | --keep ) keep=true; shift ;;
    -p | --prefix ) prefix="$2"; shift 2 ;;
    -j | --json ) json_params="$2"; shift 2 ;;
    -d | --deduplicate ) deduplicate=true; shift;;
    -t | --toi ) toi="$2"; shift 2 ;;
    -n | --dryrun ) dryrun=true; verbose=true; shift;;
    --version ) echo $VERSION ; exit 0 ;;
    -- ) shift; break ;;
    * ) break ;;
  esac
done

files=($@)
if [[ ${#files[@]} == 0 ]]; then
    usage >&2
    exit 1
fi

DEDUP="cat"
if [[ ${deduplicate} == true ]]; then
    DEDUP="dedup.py"
fi
set -o errexit  # exit on error
set -o pipefail

mkdir -p $(dirname ${prefix})

for sample in "${files[@]}"
do
    label=$(basename $sample)
    label=${label%.*}
    echocmd "${DEDUP} ${sample} >${outdir}/${label}.fa"

    echotty "Creating alignment with PRANK..." 
    echocmd "prank -d=${outdir}/${label}.fa -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"

    # callculate ancestral founders for multiple- and single-infections
    # Split the tree just below the root into two subtrees
    echocmd "treesplit.py -o  ${outdir} ${outdir}/prank.best.dnd ${outdir}/${label}.fa"

    declare -A founder
    for subtree in 'left' 'right'
    do
	# prank will not run if only given a single sequence,
	# so avoid that by testing for singletons.
	if [[ $( ${dryrun} || egrep -c '^>' ${outdir}/${subtree}.fasta) > 1 ]]
	then
	    echocmd "prank -d=${outdir}/${subtree}.fasta -o=${outdir}/${subtree} -quiet -once -f=fasta -showanc -showtree -showevents -DNA >${outdir}/prankcmd.log 2>&1"
	    founder[${subtree}]=$( ${dryrun} || echocmd prankroot.py ${outdir}/${subtree}.best.anc.dnd ${outdir}/${subtree}.best.anc.fas)
	else
	    founder[${subtree}]=$( ${dryrun} || egrep -v '^>' ${outdir}/${subtree}.fasta | tr -d '\n')
	fi
    done

    # produce the mutiple-infection ancestral sequence file
    ${dryrun} || cat <<EOF >${prefix}multiple.fa 
$(for key in "${!founder[@]}"; do printf ">%s_%s\n%s\n" ${label} ${key} ${founder[${key}]} ; done)    
EOF
    
    # single founder - process all sequences together
    # extract sequences at root of guide tree and output in FASTA format.
    ${dryrun} || cat >${prefix}single.fa <<EOF
>${label}_founder
`prankroot.py ${outdir}/prank.best.anc.dnd ${outdir}/prank.best.anc.fas`
EOF
    
	
    # inserts sample sequences into a beast configuration file.
    # parses date from the sequence ids and converts to tip date for BEAST
    echotty "Make BEAST config file..."

    # note: ${toi:+toi='${toi}'} expands to toi='${toi}' if the toi value is set, otherwise it expands to nothing.
    # http://wiki.bash-hackers.org/syntax/pe#use_an_alternate_value
    if [ -z "${json_params}" ]; then
	cmdline_params=$(for key in "${!params[@]}"; do printf "'%s'='%s' " "${key}" "${params[${key}]}"; done)
    else
	cmdline_params="--params '${json_params}'"
    fi
    echocmd "mkbeast_rv217.py --template ${TEMPLATES}/beast_strict.template ${cmdline_params}  default_rate=1.62E-2 ${toi:+toi='${toi}'} ${outdir}/prank.best.fas > ${outdir}/beast_in.xml"

    # sigh...if using the -working option of beast, the configuration file must be specified with an absolute path.
    echotty "Running BEAST..."
    echocmd "beast -working -overwrite -beagle ${outdir}/beast_in.xml"
    
    echotty 'Extracting estimated time of infection'
    echocmd "posterior_toi.py ${outdir}/beastout.log  >${prefix}toi.csv"

    # echotty "Exracting guide tree..." 
    # echocmd "treeannotator ${outdir}/beastout.trees >${outdir}/mcc.tree 2>${outdir}/treeannotator.log"

    # tr "/" "-" < ${outdir}/mcc.tree | nexus2newick.py | tr "-" "/" | tr -d "'" > ${outdir}/guidetree.tree

    # echotty "Inferring ancestral states with PRANK..." 
    # echocmd "prank -d=${outdir}/${label}.fa -t=${outdir}/guidetree.tree -o=${outdir}/prank -quiet -once -f=fasta -showanc -showtree -showevents -codon >${outdir}/prankcmd.log 2>&1"

done

