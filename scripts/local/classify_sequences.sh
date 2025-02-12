#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor
FILE=~/Desktop/datasets/test_voc/test/Alpha.fasta
CLASSIFY=$PROJECTHOME/.project_cache/trained_model.dat

cd $PROJECTHOME && julia --project Main.jl -f $FILE -o $PROJECTHOME --classify $CLASSIFY



