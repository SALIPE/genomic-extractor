#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

export JULIA_NUM_THREADS=20

PROJECT=/home/a61491/rrm-genomic-extractor/scripts/cl//verify_classes_consensus_cl.sh

PROJECT 0.002
PROJECT 0.004
PROJECT 0.006
PROJECT 0.008
PROJECT 0.01
PROJECT 0.012
PROJECT 0.014
PROJECT 0.016
PROJECT 0.018
PROJECT 0.02
PROJECT 0.03
PROJECT 0.04
PROJECT 0.05
PROJECT 0.1
PROJECT 0.15
PROJECT 0.2