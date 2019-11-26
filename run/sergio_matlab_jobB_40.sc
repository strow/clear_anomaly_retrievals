## submit 40 jobs, one per latbin (each loops over 365 timesteps, SLOW!!)
/bin/rm slurm*.out
sbatch --exclude=cnode[204,225,267] --array=1-40      --begin=now     sergio_matlab_jobB.sbatch 2
