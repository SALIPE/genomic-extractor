#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor/GREAC

# TESTDIR=~/Desktop/datasets/test_voc/test
# TRAIN=~/Desktop/datasets/test_voc/train/kmers
# GROUPNAME=covid

TESTDIR=~/Desktop/datasets/sars_cov2/test
TRAIN=~/Desktop/datasets/sars_cov2/train/kmers
GROUPNAME=covid_2

# TESTDIR=~/Desktop/datasets/bees/test
# TRAIN=~/Desktop/datasets/bees/kmers
# GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/denv/DENV_test
# TRAIN=~/Desktop/datasets/denv/DENV_train/kmers
# GROUPNAME=denv

cd $PROJECTHOME && julia --project src/GREAC.jl --group-name $GROUPNAME -w $1 benchmark --train-dir $TRAIN --test-dir $TESTDIR -m $2


