#!/bin/bash

PBS_JOBID=$(qsub scripts/run_FAST.sh)
echo "job submitted with ID $PBS_JOBID"
