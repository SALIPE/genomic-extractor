#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

export JULIA_NUM_THREADS=25

PROJECTHOME=/home/a61491/rrm-genomic-extractor
TRAIN=/home/a61491/datasets/kmers
TESTDIR=/home/a61491/datasets/test_voc/test


cd $PROJECTHOME && julia --project Main.jl classify -w $1 --test-dir $TESTDIR



