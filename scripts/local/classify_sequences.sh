#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor
FILE=~/Desktop/datasets/test_voc/train/Alpha.fasta


cd $PROJECTHOME && julia --project Main.jl -f $FILE -o teste -w 0.004 --classify



