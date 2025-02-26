#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECTHOME=/home/a61491/rrm-genomic-extractor
DIR_INPUT=/home/a61491/datasets/kmers

cd $PROJECTHOME && julia --project Main.jl extract-model -d $DIR_INPUT -w $1 



