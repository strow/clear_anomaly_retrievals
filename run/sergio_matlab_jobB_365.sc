## submit 365 jobs, one per timestep (each loops over 40 latbins, FAST!)
/bin/rm slurm*.out
sbatch --exclude=cnode[204,225,267] --array=1-30      --begin=now     sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=31-60     --begin=now+15  sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=61-90     --begin=now+30  sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=91-120    --begin=now+45  sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=121-150   --begin=now+60  sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=151-180   --begin=now+75  sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=181-210   --begin=now+90  sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=211-240   --begin=now+105 sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=241-270   --begin=now+120 sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=271-300   --begin=now+135 sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=301-330   --begin=now+150 sergio_matlab_jobB.sbatch 1
sbatch --exclude=cnode[204,225,267] --array=331-365   --begin=now+165 sergio_matlab_jobB.sbatch 1
