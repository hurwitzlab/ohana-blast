#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -p skx-normal
#SBATCH -J ohana-blast-test
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user jklynch@email.arizona.edu

module load irods

OUT_DIR="$SCRATCH/ohana-blast-test"
if [[ -d $OUT_DIR ]]; then
  rm -rf $OUT_DIR
fi

mkdir -p $OUT_DIR
iget -f /iplant/home/shared/imicrobe/test/apps/ohana-blast/test_HOT224_1_0025m.fa $OUT_DIR
ls -l $OUT_DIR

run.sh -q $OUT_DIR/test_HOT224_1_0025m.fa -o $OUT_DIR
