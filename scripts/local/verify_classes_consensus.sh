#!/bin/bash

PROJECTHOME=/home/salipe/Desktop/GitHub/rrm-genomic-extractor
DIR_INPUT=/home/salipe/Desktop/GitHub/datasets/tutorial_data/VOCs
OUTPUT=$PROJECTHOME/output_dist



cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1  -o $OUTPUT/euclidian_consensus\_$1% 



