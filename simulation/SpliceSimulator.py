import os
import argparse
import yaml
import subprocess
import shutil
import glob

def run_command(command):
    print("Running command:", " ".join(command))
    subprocess.run(command, check=True)

def load_config(config_path="config.yml"):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)
        
def parse_args():
    parser = argparse.ArgumentParser(description="Run the splicing simulation pipeline.")
    parser.add_argument("-cfig", default="config.yml",
                        help="Path to the YAML configuration file (default: config.yml)")
    return parser.parse_args()

def main():
    args = parse_args()
    basedir = os.path.dirname(os.path.abspath(__file__))
    config = load_config(args.cfig)
    os.makedirs(config['wd'], exist_ok=True)
    
    simulate_cmd = ["python", 
        f"{basedir}/src/simulate.py",
        "-incsv", config['input_csv'],
        "-ref", config['reference'],
        "-exonlabel", config['exon_label'],
        "-fastadir", f"{config['wd']}/fasta"
    ]
    run_command(simulate_cmd)
    
    polyester_cmd = ["python", 
        f"{basedir}/src/run_polyester.py",
        "-progpath", f"{basedir}/src",
        "-fastadir", f"{config['wd']}/fasta",
        "-outdir", f"{config['wd']}/simulation",
        "-lowCounts", config['lowCounts'],
        "-highCounts", config['highCounts']
    ]
    run_command(polyester_cmd)
    
    fasta2bam_cmd = ["python", 
        f"{basedir}/src/run_fasta2bam.py",
        "-indir", config['wd']
    ]
    run_command(fasta2bam_cmd)

    
    file_pattern = f"{config['wd']}/simulation*/bam/*bam*"
    destination_dir = f"{config['wd']}/bam"
    os.makedirs(destination_dir, exist_ok=True)
    source_files = glob.glob(file_pattern)
    if not source_files:
        print("No files found matching pattern.")
    for file in source_files:
        shutil.copy(file, destination_dir)
        # print(f"Copied {file} to {destination_dir}")
    

if __name__ == "__main__":
    main()