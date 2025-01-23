#!/bin/bash

PROJECTHOME=/home/salipe/Desktop/GitHub/rrm-genomic-extractor
# DIR_INPUT=/home/salipe/Desktop/GitHub/datasets/tutorial_data/VOCs
DIR_INPUT=/home/salipe/Desktop/GitHub/datasets/tutorial_data/consensus_classes
OUTPUT=$PROJECTHOME/output_dist
POSDIR=$PROJECTHOME/positions/omicron/omicron99


cd $PROJECTHOME && julia --project Main.jl -d $DIR_INPUT -w $1 -p $POSDIR -o $OUTPUT/euclidian_consensus\_$1% 



