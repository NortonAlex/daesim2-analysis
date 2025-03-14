#!/bin/bash
#PBS -N PreSegment
#PBS -l mem=186GB
#PBS -l ncpus=24
#PBS -l jobfs=24GB
#PBS -P xe2
#PBS -l walltime=10:00:00
#PBS -l storage=gdata/xe2+gdata/v10+gdata/ka08
#PBS -q normal

job_id1=$(qsub scripts/run_FAST_sensitivity_analysis_parallel.sh)
echo "job submitted with ID $job_id1"