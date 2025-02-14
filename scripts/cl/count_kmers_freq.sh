#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor
DIR_INPUT=/home/a61491/datasets/kmers
OUTPUT=$PROJECTHOME/output_freq/$1


rm -r $OUTPUT
mkdir -p $OUTPUT

cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1 -o $OUTPUT --extract-model



