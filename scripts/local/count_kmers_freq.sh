#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor
DIR_INPUT=~/Desktop/datasets/tutorial_data/kmers

cd $PROJECTHOME && julia --project Main.jl extract-model -d $DIR_INPUT -w $1 


