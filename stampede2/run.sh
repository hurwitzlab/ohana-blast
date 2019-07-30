#!/bin/bash

# Author: Ken Youens-Clark <kyclark@email.arizona.edu>

IMICROBE_WORK=/work/05066/imicrobe/iplantc.org/data

BIN=$( cd "$( dirname "$0" )" && pwd )
QUERY=""
PCT_ID=".98"
OUT_DIR="$BIN"
NUM_THREADS=$SLURM_TASKS_PER_NODE

module load blast
module load launcher
module load tacc-singularity

set -u

function lc() {
  wc -l "$1" | cut -d ' ' -f 1
}

function HELP() {
  printf "Usage:\n  %s -q QUERY -o OUT_DIR\n\n" $(basename $0)

  echo "Required arguments:"
  echo " -q QUERY"
  echo
  echo "Options:"
  echo
  echo " -p PCT_ID ($PCT_ID)"
  echo " -o OUT_DIR ($OUT_DIR)"
  echo " -n NUM_THREADS ($NUM_THREADS)"
  echo
  exit 0
}

if [[ $# -eq 0 ]]; then
  HELP
fi

while getopts :o:n:p:q:h OPT; do
  case $OPT in
    h)
      HELP
      ;;
    n)
      NUM_THREADS="$OPTARG"
      ;;
    o)
      OUT_DIR="$OPTARG"
      ;;
    p)
      PCT_ID="$OPTARG"
      ;;
    q)
      QUERY="$OPTARG"
      ;;
    :)
      echo "Error: Option -$OPTARG requires an argument."
      exit 1
      ;;
    \?)
      echo "Error: Invalid option: -${OPTARG:-""}"
      exit 1
  esac
done

if [[ $NUM_THREADS -lt 1 ]]; then
  echo "NUM_THREADS \"$NUM_THREADS\" cannot be less than 1"
  exit 1
fi

if [[ -d "$OUT_DIR" ]]; then
  mkdir -p "$OUT_DIR"
fi

BLAST_OUT_DIR="$OUT_DIR/blast-out"
if [[ ! -d "$BLAST_OUT_DIR" ]]; then
  mkdir -p "$BLAST_OUT_DIR"
fi

INPUT_FILES=$(mktemp)
if [[ -d $QUERY ]]; then
  find "$QUERY" -type f > "$INPUT_FILES"
else
  echo "$QUERY" > $INPUT_FILES
fi
NUM_INPUT=$(lc "$INPUT_FILES")

if [[ $NUM_INPUT -lt 1 ]]; then
  echo "No input files found"
  exit 1
fi

# Here is a place for fasplit.py to ensure not too
# many sequences in each query.

BLAST_DIR="$IMICROBE_WORK/ohana/blast"

if [[ ! -d "$BLAST_DIR" ]]; then
  echo "BLAST_DIR \"$BLAST_DIR\" does not exist."
  exit 1
fi

# assume this job is running on Stampede 2 SKX queue (skx-normal)
# each SKX node has 48*2=96 hardware cores
# each KNL processor has 68 cores, 4 threads per core
# TACC recommends against running 68 * 4 = 272 threads
# or we could just say 24 threads for each of 4 tasks
# see the launcher configuration below for LAUNCHER_PPN=4
BLAST_DIR="$IMICROBE_WORK/ohana/blast"
BLAST_ARGS="-outfmt 6 -num_threads 24"
BLAST_PARAM="$$.blast.param"

cat /dev/null > $BLAST_PARAM # make sure it's empty

i=0
while read INPUT_FILE; do
  BASENAME=$(basename "$INPUT_FILE")

  #dos2unix $INPUT_FILE

  let i++
  printf "%3d: %s\n" "$i" "$BASENAME"
  EXT="${BASENAME##*.}"
  TYPE="unknown"
  if [[ $EXT == 'fa'    ]] || \
     [[ $EXT == 'fna'   ]] || \
     [[ $EXT == 'fas'   ]] || \
     [[ $EXT == 'fasta' ]] || \
     [[ $EXT == 'ffn'   ]];
  then
    TYPE="dna"
  elif [[ $EXT == 'faa' ]]; then
    TYPE="prot"
  elif [[ $EXT == 'fra' ]]; then
    TYPE="rna"
  fi

  BLAST_TO_DNA=""
  if [[ $TYPE == 'dna' ]]; then
    BLAST_TO_DNA="blastn -perc_identity $PCT_ID"
  elif [[ $TYPE == 'prot' ]]; then
    BLAST_TO_DNA='tblastn'
  else
    echo "Cannot BLAST $BASENAME to DNA (not DNA or prot)"
  fi

  if [[ ${#BLAST_TO_DNA} -gt 0 ]]; then
    echo "$BLAST_TO_DNA $BLAST_ARGS -db $BLAST_DIR/contigs -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-contigs.tab" >> $BLAST_PARAM
    echo "$BLAST_TO_DNA $BLAST_ARGS -db $BLAST_DIR/genes -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-genes.tab" >> $BLAST_PARAM
  fi

  BLAST_TO_PROT=""
  if [[ $TYPE == 'dna' ]]; then
    BLAST_TO_PROT='blastx'
  elif [[ $TYPE == 'prot' ]]; then
    BLAST_TO_PROT='blastp'
  else
    echo "Cannot BLAST $BASENAME to PROT (not DNA or prot)"
  fi

  if [[ ${#BLAST_TO_PROT} -gt 0 ]]; then
    echo "$BLAST_TO_PROT $BLAST_ARGS -db $BLAST_DIR/proteins -query $INPUT_FILE -out $BLAST_OUT_DIR/$BASENAME-proteins.tab" >> $BLAST_PARAM
  fi
done < "$INPUT_FILES"
rm "$INPUT_FILES"

echo "Starting launcher for BLAST"
echo "  NUM_THREADS=$NUM_THREADS"
echo "  SLURM_JOB_NUM_NODES=$SLURM_JOB_NUM_NODES"
echo "  SLURM_NTASKS=$SLURM_NTASKS"
echo "  SLURM_JOB_CPUS_PER_NODE=$SLURM_JOB_CPUS_PER_NODE"
echo "  SLURM_TASKS_PER_NODE=$SLURM_TASKS_PER_NODE"

export LAUNCHER_DIR=$TACC_LAUNCHER_DIR
export LAUNCHER_PLUGIN_DIR=$LAUNCHER_DIR/plugins
export LAUNCHER_WORKDIR=$BIN
export LAUNCHER_RMI=SLURM
export LAUNCHER_JOB_FILE=$BLAST_PARAM

# set this one
export LAUNCHER_PPN=4

export LAUNCHER_SCHED=dynamic

echo "  LAUNCHER_PPN=$LAUNCHER_PPN"

$LAUNCHER_DIR/paramrun
echo "Ended launcher for BLAST"

rm $BLAST_PARAM

#
# Now we need to add Eggnog (and eventually Pfam, KEGG, etc.)
# annotations to the "*-genes.tab" and "*-proteins.tab" files.
#
ANNOT_PARAM="$$.annot.param"
cat /dev/null > $ANNOT_PARAM

GENE_PROTEIN_HITS=$(mktemp)
find $BLAST_OUT_DIR -size +0c -name \*-genes.tab > $GENE_PROTEIN_HITS
find $BLAST_OUT_DIR -size +0c -name \*-proteins.tab >> $GENE_PROTEIN_HITS
while read FILE; do
  BASENAME=$(basename $FILE '.tab')
  echo "Annotating $FILE"
  echo "singularity exec ohana-blast.img python3 /scripts/annotate.py -b \"$FILE\" -a \"${IMICROBE_WORK}/ohana/sqlite\" -o \"${OUT_DIR}/annotations\"" >> $ANNOT_PARAM
done < $GENE_PROTEIN_HITS

echo "Starting launcher for annotation"
export LAUNCHER_JOB_FILE=$ANNOT_PARAM

export LAUNCHER_PPN=4
# I would rather use 'dynamic' but it is failing often
export LAUNCHER_SCHED=interleaved

echo "  LAUNCHER_PPN=$LAUNCHER_PPN"

$LAUNCHER_DIR/paramrun
echo "Ended launcher for annotation"

rm $ANNOT_PARAM

#
# Now we need to extract the Ohana sequences for the BLAST hits.
#
EXTRACTSEQS_PARAM="$$.extractseqs.param"
cat /dev/null > $EXTRACTSEQS_PARAM

BLAST_HITS=$(mktemp)
find $BLAST_OUT_DIR -size +0c -name \*.tab > $BLAST_HITS
while read FILE; do
  BASENAME=$(basename $FILE '.tab')
  echo "Extracting Ohana sequences of BLAST hits for $FILE"
  echo "singularity exec ohana-blast.img python3 /scripts/extractseqs.py \"$FILE\"  \"${IMICROBE_WORK}/ohana/HOT\" \"${OUT_DIR}/ohana-hits\"" >> $EXTRACTSEQS_PARAM
done < $BLAST_HITS

echo "Starting launcher for Ohana sequence extraction"
export LAUNCHER_JOB_FILE=$EXTRACTSEQS_PARAM

export LAUNCHER_PPN=4
export LAUNCHER_SCHED=interleaved
echo "  LAUNCHER_PPN=$LAUNCHER_PPN"

$LAUNCHER_DIR/paramrun
echo "Ended launcher for Ohana sequence extraction"
rm "$EXTRACTSEQS_PARAM"

#
# Finally add a header row to the BLAST output files.
#
INSERTHDR_PARAMS="$$.inserthdr.param"
cat /dev/null > $INSERTHDR_PARAMS

BLAST_HITS=$(mktemp)
find $BLAST_OUT_DIR -size +0c -name \*.tab > $BLAST_HITS
while read FILE; do
  BASENAME=$(basename $FILE '.tab')
  echo "Inserting header in BLAST output $FILE"
  echo "singularity exec ohana-blast.img python3 /scripts/inserthdr.py \"$FILE\"" >> $INSERTHDR_PARAMS
done < $BLAST_HITS

cat "$INSERTHDR_PARAMS"

echo "Starting launcher for BLAST header insertion"
export LAUNCHER_JOB_FILE=$INSERTHDR_PARAMS

export LAUNCHER_PPN=4
export LAUNCHER_SCHED=interleaved
echo "  LAUNCHER_PPN=$LAUNCHER_PPN"

$LAUNCHER_DIR/paramrun
echo "Ended launcher for header insertion"
rm "$INSERTHDR_PARAMS"

#
# Clean up the bin directory
#
rm -rf $BIN/bin
