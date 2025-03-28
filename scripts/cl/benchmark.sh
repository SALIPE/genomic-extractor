#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor/GREAC

#TRAIN=/home/a61491/datasets/test_voc/train/kmers
#TESTDIR=/home/a61491/datasets/test_voc/test
# GROUPNAME=covid

TRAIN=/home/a61491/datasets/denv/kmers
TESTDIR=/home/a61491/datasets/denv/test
# GROUPNAME=denv



cd $PROJECTHOME && julia --project src/GREAC.jl --group-name $GROUPNAME benchmark -w $1 --train-dir $TRAIN --test-dir $TESTDIR -m $2



