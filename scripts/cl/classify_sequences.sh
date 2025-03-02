#!/bin/bash

PROJECTHOME=/home/a61491/rrm-genomic-extractor
TRAIN=/home/a61491/datasets/kmers
TESTDIR=/home/a61491/datasets/test_voc/test


cd $PROJECTHOME && julia --project Main.jl classify -w $1 --test-dir $TESTDIR



