#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor/GREAC
#TESTDIR=~/Desktop/datasets/test_voc/test
#TRAIN=~/Desktop/datasets/tutorial_data/kmers

TESTDIR=~/Desktop/datasets/denv_eq/test
TRAIN=~/Desktop/datasets/denv_eq/kmers

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache benchmark -w $1 --train-dir $TRAIN --test-dir $TESTDIR -m $2 



