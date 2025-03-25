#!/bin/bash



PROJECTHOME=~/Desktop/rrm-genomic-extractor
FILE=~/Desktop/datasets/test_voc/test

cd $PROJECTHOME && julia --project Main.jl classify --test-dir $FILE -w $1 -m kld



