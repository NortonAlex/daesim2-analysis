#!/bin/bash

PBS_JOBID=$(qsub scripts/pbs_run_FAST.sh)
echo "job submitted with ID $PBS_JOBID"
