#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor/GREAC

TESTDIR=~/Desktop/datasets/test_voc/test
TRAIN=~/Desktop/datasets/tutorial_data/kmers
GROUPNAME=covid

# TESTDIR=~/Desktop/datasets/bees/test
# TRAIN=~/Desktop/datasets/bees/kmers
# GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/denv_eq/test
# TRAIN=~/Desktop/datasets/denv_eq/kmers
# GROUPNAME=denv

cd $PROJECTHOME && julia --project src/GREAC.jl --group-name $GROUPNAME benchmark -w $1 --train-dir $TRAIN --test-dir $TESTDIR 



