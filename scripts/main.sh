#!/bin/bash

PBS_JOBID=$(qsub scripts/run_FAST_sensitivity_analysis_parallel.sh)
echo "job submitted with ID $PBS_JOBID"
