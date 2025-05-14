#!/bin/bash

PROJECTHOME=~/Desktop/rrm-genomic-extractor/GREAC

#TESTDIR=~/Desktop/datasets/test_voc/test
#TRAIN=~/Desktop/datasets/test_voc/train/kmers
#GROUPNAME=covid

# TESTDIR=~/Desktop/datasets/sars_cov2/test
# TRAIN=~/Desktop/datasets/sars_cov2/train/kmers
# GROUPNAME=covid_2

# TESTDIR=~/Desktop/datasets/bees/test
# TRAIN=~/Desktop/datasets/bees/kmers
# GROUPNAME=bees

# TESTDIR=~/Desktop/datasets/denv_eq/test
# TRAIN=~/Desktop/datasets/denv_eq/train/kmers
# GROUPNAME=denv_eq

TESTDIR=~/Desktop/rrm-genomic-extractor/comparison_scripts/castor_hiv_data/variants/test
TRAIN=~/Desktop/rrm-genomic-extractor/comparison_scripts/castor_hiv_data/variants/train/kmers
GROUPNAME=hiv

cd $PROJECTHOME && julia --project src/GREAC.jl --no-cache --group-name $GROUPNAME -w $1 benchmark --train-dir $TRAIN --test-dir $TESTDIR -m $2


