#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor
DIR_INPUT=~/Desktop/datasets/tutorial_data/kmers
OUTPUT=$PROJECTHOME/output_freq/$1


rm -r $OUTPUT
mkdir -p $OUTPUT

cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1 -o $OUTPUT --kmers-freq



