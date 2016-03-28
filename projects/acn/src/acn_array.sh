#!/bin/bash
#
#SBATCH -J sens_analysis # A single job name for the array
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue # Partition
#SBATCH --mem 1000 # Memory request
#SBATCH -t 0-00:30 # (D-HH:MM)
#SBATCH -D /n/home10/mklau/ellison_lab/fastslow/src  # Working directory
#SBATCH --mail-type=ALL # email notification
#SBATCH --mail-user=matthewklau@fas.harvard.edu

## To run this use: sbatch --array=1-5000 sens_array.sh
## Note --array=0-15%n" will limit the number of simultaneously running tasks from this job array to n.
## http://www.schedmd.com/slurmdocs/job_array.html
## https://rc.fas.harvard.edu/resources/running-jobs/

### outname="acn_mods"$(date +%Y%m%d)

outname=acn_mods20160328
Rscript acn_gmods.R ${SLURM_ARRAY_TASK_ID}
