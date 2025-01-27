#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

# source /home/a61491/.bashrc

# export JULIA_NUM_THREADS=20

PROJECTHOME=/home/a61491/rrm-genomic-extractor
# DIR_INPUT=/home/a61491/datasets/tutorial_data/VOCs
DIR_INPUT=/home/a61491/datasets/consensus
OUTPUT=$PROJECTHOME/output_dist/$1
POSDIR=$PROJECTHOME/positions/all

rm -r $OUTPUT
mkdir -p $OUTPUT

cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1 -p $POSDIR -o $OUTPUT



