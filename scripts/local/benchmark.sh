#!/bin/bash

PROJECTHOME=~/Desktop/genomic-extractor/GREAC

#TESTDIR=~/Desktop/datasets/test_voc/test
#TRAIN=~/Desktop/datasets/test_voc/train/kmers
#GROUPNAME=covid

# TESTDIR=~/Desktop/datasets/sars_cov2/test
# TRAIN=~/Desktop/datasets/sars_cov2/train/kmers
# GROUPNAME=covid_2

# TESTDIR=~/Desktop/datasets/bees/test
# TRAIN=~/Desktop/datasets/bees/kmers
# GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/dengue/test
# TRAIN=~/Desktop/datasets/dengue/train/kmers
# GROUPNAME=dengue

TESTDIR=~/Desktop/genomic-extractor/comparison_scripts/castor_hiv_data/variants/test
TRAIN=~/Desktop/genomic-extractor/comparison_scripts/castor_hiv_data/variants/train/kmers
GROUPNAME=hiv

# TESTDIR=/home/salipe/Desktop/datasets/HBV/data/test
# TRAIN=/home/salipe/Desktop/datasets/HBV/data/train/kmers
# GROUPNAME=hbv

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME -w $1 benchmark --train-dir $TRAIN --test-dir $TESTDIR -m $2


