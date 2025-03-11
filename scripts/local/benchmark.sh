#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor
#TESTDIR=~/Desktop/datasets/test_voc/test
#TRAIN=~/Desktop/datasets/tutorial_data/kmers

TESTDIR=~/Desktop/datasets/denv/kmers
TRAIN=~/Desktop/datasets/denv/kmers

cd $PROJECTHOME && julia --project Main.jl benchmark -w $1 --train-dir $TRAIN --test-dir $TESTDIR



