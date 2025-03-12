#!/bin/bash
#SBATCH -J {jName}
#SBATCH -D {wkPath}
#SBATCH -n {CoreNum}
#SBATCH -t 72:0:0
#SBATCH --mem={MEMFREE}G
#SBATCH -o {jName}.log
#SBATCH -e {jName}.log
#- End embedded arguments
echo $SLURM_JOB_NODELIST
echo 'end of HOST'
# exit
set -e

{strCMD}

echo 'DONE'