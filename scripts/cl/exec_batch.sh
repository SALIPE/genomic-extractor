#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

export JULIA_NUM_THREADS=20

PROJECT=/home/a61491/rrm-genomic-extractor/scripts/cl/benchmark.sh

$PROJECT 0.002
$PROJECT 0.004
$PROJECT 0.006
$PROJECT 0.008
$PROJECT 0.01
$PROJECT 0.05
$PROJECT 0.1
