import os

slurmScript = """#!/bin/bash
#SBATCH -J jID_jName
##SBATCH -D wkPath
#SBATCH -n CoreNum
#SBATCH --nodes=1
#SBATCH -t 72:0:0
#SBATCH --mem=MEMFREEG
#SBATCH -o wkPath/jName.o.log
#SBATCH -e wkPath/jName.e.log
#- End embedded arguments
echo $SLURM_JOB_NODELIST
module load Anaconda3
conda activate CONDAENV

RUMCMD

echo 'DONE'
"""
slurmArray = """#!/bin/bash
#SBATCH -J jID_jName
##SBATCH -D wkPath
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=CoreNum
#SBATCH -t 72:0:0
#SBATCH --mem=MEMFREEG
#SBATCH --array=1-ARRAYN
#SBATCH -o wkPath/jName.o.log
#SBATCH -e wkPath/jName.e.log
#- End embedded arguments
echo $SLURM_JOB_NODELIST
module load Anaconda3
conda activate CONDAENV

#echo "Task: $SLURM_ARRAY_TASK_ID"
RUMCMD

echo "DONE: $SLURM_ARRAY_TASK_ID"
"""


def create_sbatch_script(jName, output_path, cmd, n_tasks=4, mem_free=None, time_limit="72:0:0"):
    if mem_free is None:
        mem_free = 16 * n_tasks
    mem_str = f"{mem_free}G"
    log_file = f"{jName}.log"
    cmd_str = "\n".join(cmd)
    script = f"""#!/bin/bash
    #SBATCH -J {jName}
    #SBATCH -D {output_path}
    #SBATCH -n {n_tasks}
    #SBATCH -t {time_limit}
    #SBATCH --mem={mem_str}
    #SBATCH -o {log_file}
    #SBATCH -e {log_file}
    #- End embedded arguments
    echo $SLURM_JOB_NODELIST
    module load Anaconda3
    conda activate CONDAENV

    {cmd_str}

    echo 'DONE'
    """
    return script

def write_sbatch_script(script, output_path, jName):
    script_path = os.path.join(output_path, jName + '.sh')
    with open(script_path, 'w') as f:
        f.write(script)
    return script_path