#!/bin/bash

#$ -o /home/a61491/.outputs
#$ -e /home/a61491/.errs

source /home/a61491/.bashrc

PROJECT=/home/a61491/rrm-genomic-extractor/scripts/cl/benchmark.sh

$PROJECT 0.002
$PROJECT 0.004
$PROJECT 0.006
$PROJECT 0.008
$PROJECT 0.01
