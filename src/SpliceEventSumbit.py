import os
import subprocess
import argparse
import yaml

class SpliceHarmJobSubmitter:
    """
    A class to handle the creation and submission of Slurm jobs for splicing harmonization.
    Configuration parameters such as working directory, reference file, and data directory
    are loaded from a YAML file.
    """
    def __init__(self, config_path):
        with open(config_path, "r") as f:
            config = yaml.safe_load(f)
        
        # Load configuration parameters
        self.WD = config.get("wdir")
        self.reference = config.get("genome_gtf")
        self.data_dir = config.get("indir")
        self.create_script_path = os.path.join(os.path.dirname(os.path.dirname(os.path.realpath(__file__))),  'src', 'SpliceEvent', 'create_scripts.py')
        self.majiq_cutoff_val = config.get("majiq_cutoff_val", 0.95)
        self.junction_filter = config.get("junction_filter", False)
        
        if not (self.WD and self.reference and self.data_dir and self.create_script_path):
            raise ValueError("Configuration file is missing required parameters.")

    def create_job_scripts(self):
        create_scripts_cmd = [
            "python3",
            self.create_script_path,
            "-indir", self.data_dir,
            "-wdir", self.WD,
            "-ref", self.reference,
            "-majiq_cutoff_val", str(self.majiq_cutoff_val),
            "-junction_filter", str(self.junction_filter)
        ]
        try:
            subprocess.run(create_scripts_cmd, check=True)
            print("Job scripts created successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error creating job scripts: {e}")
            raise

    def submit_job(self, script_name, dependency_id=None):
        """
        Submits a job script to Slurm via sbatch. If a dependency job ID is provided,
        the script is submitted with a dependency on that job's completion.
        
        Returns:
            job_id (str): The job ID returned from sbatch.
        """
        cmd = ["sbatch"]
        if dependency_id:
            cmd.append(f"--dependency=afterok:{dependency_id}")
        cmd.append(script_name)
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            output = result.stdout.strip()
            # Expected output format: "Submitted batch job <job_id>"
            job_id = output.split()[-1]
            print(f"Submitted {script_name} with job ID: {job_id}")
            return job_id
        except subprocess.CalledProcessError as e:
            print(f"Error submitting {script_name}: {e.stderr}")
            return None

    def submit_jobs(self):
        """
        """
        os.chdir(self.WD)
        subfolders = sorted([f for f in os.listdir(self.WD) if os.path.isdir(os.path.join(self.WD, f))])
        print(f"Found {len(subfolders)} subdirectories.")
        
        for folder in subfolders:
            folder_path = os.path.join(self.WD, folder)
            run1_path = os.path.join(folder_path, "run1.sh")
            
            if os.path.isfile(run1_path):
                print(f"Processing subdirectory: {folder}")
                os.chdir(folder_path)
                
                # Submit run1.sh
                job_id1 = self.submit_job("run1.sh")
                
                if job_id1 and os.path.isfile("run2.sh"):
                    # Submit run2.sh with dependency on run1.sh
                    job_id2 = self.submit_job("run2.sh", dependency_id=job_id1)
                    
                    if job_id2 and os.path.isfile("run3.sh"):
                        # Submit run3.sh with dependency on run2.sh
                        job_id3 = self.submit_job("run3.sh", dependency_id=job_id2)
                        
                        if job_id3 and os.path.isfile("run4.sh"):
                            # Submit run4.sh with dependency on run3.sh
                            job_id4 = self.submit_job("run4.sh", dependency_id=job_id3)
                            
                            if job_id4 and os.path.isfile("run5.sh"):
                                # Submit run5.sh with dependency on run4.sh
                                self.submit_job("run5.sh", dependency_id=job_id4)
                # Return to the main working directory for the next subfolder.
                os.chdir(self.WD)
            else:
                print(f"No 'run1.sh' found in {folder}")

    def run(self):
        """
        Runs the entire job submission process: creating job scripts and submitting them.
        """
        try:
            self.create_job_scripts()
        except subprocess.CalledProcessError:
            print("Failed to create job scripts. Exiting.")
            return
        self.submit_jobs()

    @classmethod
    def main(cls, config_path):
        """
        The class method to initiate the job submission process.
        Usage:
            SpliceHarmJobSubmitter.main(config_path)
        """
        job_submitter = cls(config_path)
        job_submitter.run()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SpliceHarm Job Submitter")
    parser.add_argument("-cfig", required=True, help="Path to the config.yml file")
    args = parser.parse_args()
    
    # Call the class method main with the provided configuration file path.
    SpliceHarmJobSubmitter.main(args.cfig)
