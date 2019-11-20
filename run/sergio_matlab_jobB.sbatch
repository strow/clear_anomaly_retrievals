#!/bin/bash

# run this with   sbatch --constraint=hpcf2013       --array=1-365 sergio_matlab_jobB.sbatch  for clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomaly
#   or            sbatch --exclude=cnode[204,225,260,267] --array=1-365 sergio_matlab_jobB.sbatch 1
#
# run this with   sbatch --constraint=hpcf2013       --array=1-40 sergio_matlab_jobB.sbatch  for clust_run_retrieval_latbins_AIRS_iasitimespan_loop_anomalyV2
#   or            sbatch --exclude=cnode[204,225,260,267] --array=1-40 sergio_matlab_jobB.sbatch 2
#
# check stats after running using eg
# sacct --name=KCARTA_BREAKOUTS --state=CD -S 2016-05-07 --format=elapsed,maxrss
#
#  Name of the job:
#SBATCH --job-name=ANOMALY

#  N specifies that 1 job step is to be allocated per instance of matlab
#SBATCH -N1

#  This specifies the number of cores per matlab session will be
#available for parallel jobs
#SBATCH --cpus-per-task 1

#  Specify the desired partition develop/batch/prod
#SBATCH --partition=batch
##SBATCH --partition=high_mem

#  Specify the qos and run time (format:  dd-hh:mm:ss)
##SBATCH --qos=medium+
#SBATCH --qos=short+

#SBATCH --time=0:59:00 

##  This is in MB, less aggressive
##SBATCH --mem-per-cpu=12000
##  This is in MB, very lean, don;t go below this as java starts crying
##SBATCH --mem-per-cpu=4000
##  This is in MB, very aggressive but I have been running outta memory, CRIS convolve
##SBATCH --mem-per-cpu=24000

##  This is in MB, easy peasy, 18000 works
##SBATCH --mem-per-cpu=18000
##SBATCH --mem-per-cpu=8000
##SBATCH --mem-per-cpu=6000 
#SBATCH --mem-per-cpu=6000 

##  Specify the job array (format:  start-stop:step) 
# AIRS
#srun matlab -nojvm -nodisplay -r "clust_convert_strowrates2oemrates_anomaly; exit"
#srun matlab -nojvm -nodisplay -nosplash -r "clust_run_retrieval_latbins_AIRS_loop_anomaly; exit"
#
# IASI
#srun matlab -nojvm -nodisplay -r "clust_convert_strowrates2oemrates_anomaly_11year_iasitimespan; exit"
#

if [ $# -gt 0 ]; then
  echo "Your command line contains $# arguments"
elif [ $# -eq 0 ]; then
  echo "Your command line contains no arguments"
fi

if [[ "$1" -eq "" ]]; then
  # this is fixed timestep, loop over 40 latbins 
  echo "cmd line arg = DNE, clust_run_retrieval_latbinsloop_anomaly => JOB = timestep"
  srun matlab -nodisplay -nosplash -r "addpath ../src/; clust_run_retrieval_latbins_loop_anomaly; exit"
elif [[ "$1" -eq "1" ]]; then
  # this is fixed timestep, loop over 40 latbins 
  echo "cmd line arg = 1, clust_run_retrieval_latbinsloop_anomaly => JOB = timestep"
  srun matlab -nodisplay -nosplash -r "addpath ../src/; clust_run_retrieval_latbins_loop_anomaly; exit"
elif [[ "$1" -eq "2" ]]; then
  # this fixed latbin, loop over 365 timesteps
  echo "cmd line arg = 2, clust_run_retrieval_latbinsloop_anomalyV2 => JOB = latbin"
  srun matlab -nodisplay -nosplash -r "addpath ../src/; clust_run_retrieval_latbins_loop_anomalyV2; exit"
fi
