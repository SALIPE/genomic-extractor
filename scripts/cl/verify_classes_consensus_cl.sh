#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

export JULIA_NUM_THREADS=16

PROJECTHOME=/home/a61491/rrm-genomic-extractor
DIR_INPUT=/home/a61491/datasets/tutorial_data/VOCs
OUTPUT=$PROJECTHOME/output_convergence

cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1  -o $OUTPUT/euclidian_consensus\_$1% 



