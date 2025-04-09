#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor/GREAC
DATAHOME=/tmp2/felipe

# TRAIN=/$HOME/test_voc/train/kmers
# TESTDIR=/$HOME/test_voc/test
# GROUPNAME=covid

TRAIN=$DATAHOME/sars_cov2/train/kmers
TESTDIR=$DATAHOME/sars_cov2/test
GROUPNAME=covid_2

# TESTDIR=$HOME/bees/test
# TRAIN=$HOME/bees/kmers
# GROUPNAME=bees

# TRAIN=/$HOME/denv/kmers
# TESTDIR=/$HOME/denv/test
# GROUPNAME=denv

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME -w $1 benchmark  --train-dir $TRAIN --test-dir $TESTDIR -m $2

