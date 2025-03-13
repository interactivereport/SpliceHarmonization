import os
import subprocess
from subprocess import call
import argparse
import glob
from pathlib import Path

def fasta_to_bam(fasta_dir, output_dir, python_script, cond_prefix, base_name):
    """
    Converts all FASTA files in a directory to FASTQ files using an external Python script.

    Parameters:
    fasta_dir (str): The directory containing FASTA files.
    output_dir (str): The directory where FASTQ files will be stored.
    python_script (str): Path to the Python script that performs the conversion.
    """
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    qc_dir = os.path.join(os.path.dirname(fasta_dir), 'qc')
    bam_dir = os.path.join(os.path.dirname(fasta_dir), 'bam')
    os.makedirs(qc_dir, exist_ok=True)
    # os.makedirs(bam_dir, exist_ok=False)

    # Iterate over each FASTA file in the directory
    for filename in os.listdir(fasta_dir):
        if filename.endswith(".fasta"):
            sample_name = os.path.splitext(filename)[0]
            fasta_path = os.path.join(fasta_dir, filename)
            fastq_path = os.path.join(output_dir, f"{cond_prefix}_{sample_name}.fastq")

            # Call the external Python script for conversion
            call(["python", python_script, fasta_path, fastq_path])

    fastq_files = glob.glob(output_dir + '/*.fastq')

    paired_files = {}
    for file in fastq_files:
        sample_key = '_'.join(file.split('/')[-1].split('_')[:3])
        if sample_key in paired_files:
            paired_files[sample_key].append(file)
        else:
            paired_files[sample_key] = [file]

    for key in paired_files.keys():
        paired_files[key].sort()

    
    commands = []
    for pair in paired_files.values():
        if len(pair) == 2:  # Ensure there are two files in the pair
            cmd = f"/bin/bash -l -c 'module load FastQC/0.12.1-Java-11 && fastqc {' '.join(str(p) for p in pair)} --outdir {qc_dir} -t 4'"
            commands.append(cmd)
    for cmd in commands:
        try:
            # Using subprocess.run to execute the command
            result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
            print(f"Command executed successfully: {cmd}")
            print(f"Output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            # Error handling
            print(f"Failed to execute command: {cmd}")
            print(f"Error: {e.stderr}")


    fastqc_cmd = f"/bin/bash -l -c 'module load Anaconda3/2022.05 && conda activate /edgehpc/apps/gb/anaconda3/4.9.2/envs/multiqc && cd {qc_dir} && multiqc -n multiqc_report_fastq {qc_dir}'"

    try:
        result = subprocess.run(fastqc_cmd, shell=True, check=True, capture_output=True, text=True)
        print("MultiQC executed successfully:")
        print("Output:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Failed to execute MultiQC:")
        print("Error:", e.stderr)



    # ref_sampdir = "/mnt/depts/dept04/compbio/projects/TST12320/4simulation_2/simulated_data/CNTL_fastq"
    # star_dir = "/mnt/depts/dept04/compbio/projects/TST12320/4simulation_2/simulated_data/CNTL_fastq/bam"
    # star_seq_fastq_samples = ["CNTL_sample_01", "CNTL_sample_02", "CNTL_sample_03"]  # Example sample names
    genomeDir='/home/ychen12/splicing/TST12320/4simulation_2/STAR/STARindex'
    # Loop over samples to process each one
    for i, sample in enumerate(paired_files.keys()):
    
        out_prefix = Path(bam_dir, sample, sample)
        Path(bam_dir, sample).mkdir(parents=True, exist_ok=True)
        fastq_1 = Path(output_dir, f"{sample}_1.fastq")  
        fastq_2 = Path(output_dir, f"{sample}_2.fastq")  

        star_seq_cmd = ("/bin/bash -l -c '"
            "module use /usr/prog/modules/all; "
            "module purge; "
            "module load STAR/2.7.2d-GCC-12.2.0; "
            "STAR --runThreadN 8 "
            "--runMode alignReads "
            f"--genomeDir {genomeDir} "
            f"--readFilesIn {fastq_1} {fastq_2} "
            f"--outFileNamePrefix {out_prefix} "
            "--outSAMtype BAM Unsorted "
            "--outReadsUnmapped Fastx'"
        )
        
        
        result = subprocess.run(star_seq_cmd, shell=True, capture_output=True, text=True)

        # Optionally, print the output or handle errors
        if result.returncode == 0:
            print(f"STAR alignment successful for sample {sample}")
        else:
            print(f"Error in STAR alignment for sample {sample}: {result.stderr}")


        in_bam = f"{bam_dir}/{sample}/{sample}Aligned.out.bam"
        out_bam = f"{bam_dir}/{cond_prefix}_{base_name}_{i}.sorted_Aligned.out.bam"
        star_seq_cmd = ("/bin/bash -l -c '"
            "module use /usr/prog/modules/all; "
            "module purge; "
            "module load SAMtools/1.9-GCC-12.2.0; "
            f"samtools sort {in_bam} -o {out_bam}; "
            f"samtools index {out_bam}; "
            "'" 
        )
    
        result = subprocess.run(star_seq_cmd, shell=True, capture_output=True, text=True)

        if result.returncode == 0:
            print(f"STAR alignment successful for sample {sample}")
        else:
            print(f"Error in STAR alignment for sample {sample}: {result.stderr}")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Check for splice event annotation in CSV.")
    # parser.add_argument('-program', dest="program", help ='program directiony', required=True)
    parser.add_argument('-indir', dest='indir', help="input directionary", required=True)
    args = parser.parse_args()
    python_script = "/home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/fasta_to_fastq.py"
    base_folder_paths = glob.glob(args.indir + '/simulation*/')


    for base_dir in base_folder_paths:
        normalized_path = base_dir.rstrip('/')
        base_name = os.path.basename(normalized_path)

        fasta_dir = os.path.join(base_dir, 'condition_alt')
        output_dir = os.path.join(base_dir, 'ALT_fastq')
        fasta_to_bam(fasta_dir, output_dir, python_script, 'ALT', base_name)
        fasta_dir = os.path.join(base_dir, "condition_cntl")
        output_dir = os.path.join(base_dir, "CNTL_fastq")
        fasta_to_bam(fasta_dir, output_dir, python_script, 'CNTL', base_name)