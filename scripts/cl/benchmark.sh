#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor/GREAC
HOME=/home/a61491/datasets

# TRAIN=/$HOME/test_voc/train/kmers
# TESTDIR=/$HOME/test_voc/test
# GROUPNAME=covid

TESTDIR=$HOME/bees/test
TRAIN=$HOME/bees/kmers
GROUPNAME=bees

# TRAIN=/$HOME/denv/kmers
# TESTDIR=/$HOME/denv/test
# GROUPNAME=denv



cd $PROJECTHOME && julia --project src/GREAC.jl --group-name $GROUPNAME benchmark -w $1 --train-dir $TRAIN --test-dir $TESTDIR #-m $2



