#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor
TRAIN=/home/a61491/datasets/test_voc/train/kmers
TESTDIR=/home/a61491/datasets/test_voc/test

#TRAIN=/home/a61491/datasets/denv/kmers
#TESTDIR=/home/a61491/datasets/denv/DENV_test


cd $PROJECTHOME && julia --project Main.jl benchmark -w $1 --train-dir $TRAIN --test-dir $TESTDIR



