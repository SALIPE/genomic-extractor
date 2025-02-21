#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor
FILE=~/Desktop/datasets/test_voc/test
CLASSIFY=$PROJECTHOME/.project_cache/extracted_features.dat

cd $PROJECTHOME && julia --project Main.jl classify -f $FILE  



