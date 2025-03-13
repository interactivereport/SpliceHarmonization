import argparse,os,shutil,glob,ast,time,configparser,yaml,pickle,re,time,sys,subprocess
import pandas as pd
from datetime import datetime
from io import StringIO
from tabulate import tabulate
from slurmUtils import slurmScript, slurmArray

#!/usr/bin/env python
import argparse
import os, shutil, glob, time, configparser, yaml, pickle, re, sys, subprocess
import pandas as pd
from datetime import datetime
from io import StringIO
from tabulate import tabulate

class SpliceEventPipeline:
    def __init__(self, config):
        """
        Initialize the pipeline using settings loaded from the YAML configuration.
        Expected keys include:
          - indir: Input data directory
          - wdir: Working directory (output)
          - ref: Reference GTF (optional; if not provided, it is deduced)
          - bamindir: Directory with bam files (if junction_filter is used)
          - specific_comparison: Comma-separated string to filter comparisons
          - junction_filter: Boolean flag to enable junction filtering
          - majiq_cutoff_val: Cutoff value for majiq filtering (default 0.95)
          - junctionFC: Fold-change threshold (default 1.2)
          - junctionMAX: Maximum count (default 25)
          - samplesheet: Path to sample sheet file
          - addexons: Path to added exons sheet file (optional)
          - cmdonly: If True, only generate SLURM scripts without submission
          - retryN: Number of times to retry a job (default 2)
          - GTFpath: Base path to GTF files (default provided)
        """
        self.config = config
        self.indir = config.get("indir")
        self.wdir = config.get("wdir")
        self.ref = config.get("genome_gtf")
        self.bamindir = config.get("bamindir")
        self.specific_comparison = config.get("specific_comparison")
        self.junction_filter = config.get("junction_filter", False)
        self.majiq_cutoff_val = config.get("majiq_cutoff_val", 0.95)
        self.junctionFC = config.get("junctionFC", 1.2)
        self.junctionMAX = config.get("junctionMAX", 25)
        self.samplesheet = config.get("samplesheet")
        self.addexons = config.get("addexons", os.path.join(self.get_pipe_path(), 'data', 'STMN2_add_exon.csv'))
        self.cmdonly = config.get("cmdonly", False)
        self.retryN = config.get("retryN", 2)
        self.GTFpath = config.get("GTFpath", '/edgehpc/dept/compbio/reference/DNAnexus_references')
        
        # Determine the pipeline base path from the current file location.
        self.strPipePath = self.get_pipe_path()
        self.validate_directories()
        
        # If no reference is provided, try to deduce it from an analysis file.
        if self.ref is None:
            analysis_files = glob.glob(os.path.join(self.indir, 'analysis-*'))
            if analysis_files:
                strAna = analysis_files[0]
                if os.path.isfile(strAna):
                    with open(strAna, 'r') as f:
                        aInfo = yaml.safe_load(f)
                    version = aInfo["FASTR_REFERENCE_VERSION"]
                    species = aInfo["FASTR_REFERENCE_SPECIES"]
                    self.ref = os.path.join(self.GTFpath, 'rnaseq', species, version, f"{version}.transcript.gtf.gz")
        if self.ref is None or not os.path.isfile(self.ref):
            self.exit_with_msg(f"Missing reference GTF: {self.ref}")
        
        # SLURM script templates
        self.slurmScript = slurmScript
        self.slurmArray = slurmArray
    
    def get_pipe_path(self):
        """Return the base path of the pipeline."""
        return os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    
    def exit_with_msg(self, msg=""):
        if len(msg) > 3:
            print(msg)
        self.print_msg_power()
        sys.exit(1)
    
    def print_msg_power(self):
        print("\nPowered by the Research Data Sciences Group [yirui.chen@biogen.com;zhengyu.ouyang@biogen.com]")
        print("------------")
    
    def print_msg_init(self):
        print("\n\n*****", datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "*****")
        git_dir = os.path.join(self.strPipePath, ".git")
        if os.path.isdir(git_dir):
            gitConfig = configparser.ConfigParser()
            gitConfig.read(os.path.join(git_dir, "config"))
            url = gitConfig.get('remote "origin"', 'url', fallback="Unknown URL")
            git_log_path = os.path.join(git_dir, "logs", "HEAD")
            gitLog = pd.read_csv(git_log_path, sep="\t", header=None)
            parts = gitLog.iloc[-1, 0].split(" ")
            info_lines = [
                f"###########\n## SpliceEvent: {url}\n",
                f"## Run Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
                f"## Pipeline Path: {self.strPipePath}\n",
                f"## Pipeline Date: {datetime.fromtimestamp(int(parts[-2])).strftime('%Y-%m-%d %H:%M:%S')} {parts[-1]}\n",
                f"## git HEAD: {parts[1]}\n###########\n\n",
                f"{self.strPipePath}/SpliceEvent_run {' '.join(sys.argv[1:])}\n"
            ]
            info_str = ''.join(info_lines)
            print(info_str)
            with open(os.path.join(self.wdir, '_cmdline_' + datetime.now().strftime('%Y%m%d.%H%M%S')), 'w') as f:
                f.write(info_str)
    
    def validate_directories(self):
        if not os.path.isdir(self.indir):
            self.exit_with_msg(f"Missing fastr splicing folder: {self.indir}")
        if not os.path.isdir(os.path.dirname(self.wdir)):
            self.exit_with_msg(f"Missing output folder: {os.path.dirname(self.wdir)}")
        os.makedirs(self.wdir, exist_ok=True)
    
    def get_comparison(self):
        stringtie_path = os.path.join(self.indir, "stringtie")
        if not os.path.isdir(stringtie_path):
            self.exit_with_msg(f"Missing stringtie folder: {stringtie_path}")
        comps = [re.sub('.combined.gtf$', '', os.path.basename(x))
                 for x in glob.glob(os.path.join(stringtie_path, "*.combined.gtf"))]
        if self.specific_comparison is not None:
            speComp = self.specific_comparison.split(',')
            comps = [comp for comp in comps if any(s in comp for s in speComp)]
        if len(comps) == 0:
            self.exit_with_msg("No specific comparison found from stringtie")
        return comps
    
    def get_one_cmd(self, oneComp):
        """
        Build a dictionary of SLURM commands for one comparison.
        """
        data_folder = self.indir
        bam_folder = self.bamindir
        wd_folder = self.wdir
        reference = self.ref
        sample_sheet_file = self.samplesheet
        file_prefix = os.path.join(wd_folder, oneComp)
        os.makedirs(os.path.join(file_prefix, 'junction_prep'), exist_ok=True)

        comparisonRef_name = oneComp.split('_vs_')[-1]
        rmats_input = os.path.join(data_folder, 'rmats', oneComp + '_*JCEC.txt')
        leafcutter_input = os.path.join(data_folder, 'leafcutter', oneComp + '*_cluster_significance.txt')
        majiq_input = os.path.join(data_folder, 'majiq', oneComp + '*voila.tsv')
        rmats_output = os.path.join(file_prefix, 'junction_prep', 'rmats_junction_prep.csv')
        leafcutter_output = os.path.join(file_prefix, 'junction_prep', 'leafcutter_junction_prep.csv')
        majiq_output = os.path.join(file_prefix, 'junction_prep', 'majiq_junction_prep.csv')
        junction_input = os.path.join(file_prefix, 'junction_prep', '{method}_junction_prep.csv')
        sj_output = os.path.join(file_prefix, 'junction_prep', 'sj.csv')
        annotation = os.path.join(data_folder, 'stringtie', oneComp + '.combined.gtf')
        harm_prep_dir = os.path.join(file_prefix, 'harm')
        mxe_prep_dir = os.path.join(file_prefix, 'harm', 'mxe')
        harm_post_input = os.path.join(harm_prep_dir, 'harm*.csv')
        mxe_post_input = os.path.join(mxe_prep_dir, 'harm*.csv')
        harm_post_dir = os.path.join(file_prefix, 'harm', 'out')
        mxe_post_dir = os.path.join(file_prefix, 'harm', 'mxe', 'out')
        relabel = os.path.join(harm_prep_dir, 'exon_relabel.csv')
        majiq_cutoff = os.path.join(harm_prep_dir, 'majiq_cutoff_value.csv')
        final_output = os.path.join(file_prefix, 'out')
        final_mxe_output = os.path.join(file_prefix, 'out', 'mxe')
        final_allgene_output = os.path.join(file_prefix, 'out', 'all_gene')
        all_gene_input = os.path.join(final_allgene_output, '*.csv')
        all_junction_output = os.path.join(file_prefix, 'out', 'all_gene_junction')
        
        os.makedirs(harm_post_dir, exist_ok=True)
        os.makedirs(mxe_post_dir, exist_ok=True)
        os.makedirs(final_mxe_output, exist_ok=True)
        os.makedirs(final_allgene_output, exist_ok=True)
        os.makedirs(all_junction_output, exist_ok=True)
        
        bam_pattern = "Aligned.sortedByCoord.out.bam"
        altbam_files, refbam_files = None, None
        if self.junction_filter and sample_sheet_file:
            sample_sheet = pd.read_csv(sample_sheet_file, sep="\t")
            sample_sheet["SampleName"] = sample_sheet["SampleName"].astype(str)
            sample_sheet["Groups"] = sample_sheet["Groups"].astype(str)
            alt_group = file_prefix.split('_vs_')[0]
            ref_group = file_prefix.split('_vs_')[1]
            alt_rows = sample_sheet[sample_sheet['Groups'] == alt_group].index
            print(alt_group)
            alt_bams = sample_sheet.SampleName[alt_rows] + '.' + bam_pattern
            print(alt_bams)
            ref_rows = sample_sheet[sample_sheet['Groups'] == ref_group].index
            ref_bams = sample_sheet.SampleName[ref_rows] + '.' + bam_pattern
            if self.junction_filter and bam_folder:
                altbam_files = [os.path.join(bam_folder, x) for x in alt_bams.to_list()]
                refbam_files = [os.path.join(bam_folder, x) for x in ref_bams.to_list()]
        inputjson = os.path.join(final_output, 'libdepth_counts.json')
        
        allCMD = {}
        # run 1
        cmds_run1 = []
        cmds_run1.append(f"Rscript {os.path.join(self.strPipePath, 'src', 'rmats_junction_prep.R')} -i '{rmats_input}' -o {rmats_output} -n {comparisonRef_name}")
        cmds_run1.append(f"Rscript {os.path.join(self.strPipePath, 'src', 'leafcutter_junction_prep.R')} -i '{leafcutter_input}' -o {leafcutter_output}")
        cmds_run1.append(f"Rscript {os.path.join(self.strPipePath, 'src', 'majiq_junction_prep.R')} -i '{majiq_input}' -o {majiq_output}")
        cmds_run1.append(f"Rscript {os.path.join(self.strPipePath, 'src', 'annotation_SJ.R')} -a {annotation} -r {reference} -i {junction_input} -o {sj_output}")
        extra = f" -addexons {self.addexons}" if self.addexons else ""
        cmds_run1.append(f"python -u {os.path.join(self.strPipePath, 'src', 'harm_prep_remove_mapping_add_exons.py')} -a {annotation} -r {reference} -sj {sj_output} -o {harm_prep_dir}{extra}")
        allCMD['run1'] = {'cmd': '\n'.join(cmds_run1), 'wd': file_prefix, 'mem': 96, 'cpu': 8, 'task': 1, 'run': 0}
        
        # run 2
        cmds_run2 = []
        cmds_run2.append(f"combined_files=($(ls {harm_post_input}) $(ls {mxe_post_input}))")
        cmds_run2.append("input_csv=${combined_files[$SLURM_ARRAY_TASK_ID-1]}")
        cmds_run2.append(f"python -u {os.path.join(self.strPipePath, 'src', 'harm_assign_job_array.py')} -exonlabel {relabel} -infile $input_csv -outdir {harm_post_dir} -majiqcutofffile {majiq_cutoff}")
        allCMD['run2'] = {'cmd': '\n'.join(cmds_run2), 'wd': file_prefix, 'mem': 16, 'cpu': 1, 'task': 80, 'run': 0}
        
        # run 3
        cmds_run3 = []
        cmds_run3.append(f"python -u {os.path.join(self.strPipePath, 'src', 'post_process.py')} -indir {harm_post_dir} -outdir {final_output} -allgene {final_allgene_output} -majiq_score_cutoff {self.majiq_cutoff_val}")
        cmds_run3.append(f"python -u {os.path.join(self.strPipePath, 'src', 'post_process_4mxe.py')} -indir {mxe_post_dir} -outdir {final_mxe_output}")
        if self.junction_filter:
            cmds_run3.append(f"python -u {os.path.join(self.strPipePath, 'src', 'libdepth_count.py')} -outdir {final_output} -Altbam {altbam_files} -Refbam {refbam_files}")
        allCMD['run3'] = {'cmd': '\n'.join(cmds_run3), 'wd': file_prefix, 'mem': 96, 'cpu': 8, 'task': 1, 'run': 0}
        
        # run 4 (and run 5 if junction_filter)
        cmds_run4 = []
        if self.junction_filter:
            cmds_run4.append(f"combined_files=($(ls {all_gene_input}) {os.path.join(final_mxe_output, 'all_cleaned.csv')})")
            cmds_run4.append("files_per_task=$(( (${#combined_files[@]} + SLURM_ARRAY_TASK_MAX - 1) / SLURM_ARRAY_TASK_MAX ))")
            cmds_run4.append("start_idx=$(( (SLURM_ARRAY_TASK_ID - 1) * files_per_task ))")
            cmds_run4.append("end_idx=$(( start_idx + files_per_task - 1 ))")
            cmds_run4.append("if [ $end_idx -ge ${#combined_files[@]} ]; then\n\tend_idx=$((${#combined_files[@]} - 1))\nfi")
            cmds_run4.append('for i in $(seq $start_idx $end_idx); do\n\tinput_csv="${combined_files[$i]}"')
            cmds_run4.append(f'\tpython -u {os.path.join(self.strPipePath, "src", "post_process_MXE_junction.py")} -infile $input_csv -outdir {all_junction_output} -Altbam "{altbam_files}" -Refbam "{refbam_files}" -inputjson {inputjson} -FCthreshold {self.junctionFC} -maxcount {self.junctionMAX}\ndone')
            allCMD['run4'] = {'cmd': '\n'.join(cmds_run4), 'wd': file_prefix, 'mem': 8, 'cpu': 1, 'task': 50, 'run': 0}
            cmds_run5 = []
            cmds_run5.append(f"python -u {os.path.join(self.strPipePath, 'src', 'post_process_MXE_junction.py')} -indir {all_junction_output} -mxeindir {final_mxe_output} -outdir {final_mxe_output}")
            allCMD['run5'] = {'cmd': '\n'.join(cmds_run5), 'wd': file_prefix, 'mem': 16, 'cpu': 1, 'task': 1, 'run': 0}
        else:
            cmds_run4.append(f"python -u {os.path.join(self.strPipePath, 'src', 'post_process_MXE.py')} -indir {final_allgene_output} -mxeindir {final_mxe_output} -outdir {final_mxe_output}")
            allCMD['run4'] = {'cmd': '\n'.join(cmds_run4), 'wd': file_prefix, 'mem': 96, 'cpu': 8, 'task': 1, 'run': 0}
        
        return allCMD
    
    def parallel_monitor(self, cmdList):
        stat = {}
        # First pass: check tasks
        for k in cmdList:
            for k1 in sorted(cmdList[k].keys()):
                if cmdList[k][k1]['run'] < 0:
                    continue
                if self.check_task(cmdList[k][k1], k1):
                    break
        # Build status summary
        stat = {}
        running = False
        for k in cmdList:
            wait = False
            for k1 in sorted(cmdList[k].keys()):
                if k1 not in stat:
                    stat[k1] = {}
                run_val = cmdList[k][k1]['run']
                if wait:
                    stat[k1][k] = "Waiting"
                    continue
                if run_val == -2:
                    stat[k1][k] = "Resumed"
                elif run_val == -1:
                    stat[k1][k] = "Completed"
                elif run_val > self.retryN:
                    stat[k1][k] = "Failed"
                else:
                    stat[k1][k] = f"Running_{cmdList[k][k1].get('jID', '')}_{run_val}"
                    running = True
                if stat[k1][k].startswith('Running') or stat[k1][k].startswith('Failed'):
                    wait = True
        with open(os.path.join(self.wdir, 'status.log'), 'w') as f:
            f.write(tabulate(pd.DataFrame(stat), headers='keys', tablefmt='plain') + "\n\n")
        return running
    
    def check_task(self, aTask, rID):
        aID = []
        if 'jID' in aTask:
            jStat, aID = self.get_slurm_job_status(aTask['jID'])
            if jStat == "COMPLETED":
                aTask['run'] = -1
                return False
            elif jStat == "RUNNING":
                return True
            elif jStat == "FAILED":
                e_log = os.path.join(aTask['wd'], rID + ".e.log")
                if os.path.isfile(e_log):
                    shutil.move(e_log, os.path.join(aTask['wd'], f"{rID}_{aTask['run']}.e.log"))
                with open(os.path.join(aTask['wd'], rID + ".o.log"), 'a') as f:
                    f.write("\nFailed with at least some jobs!")
                aTask['run'] += 1
            elif jStat == 'CANCELLED':
                with open(os.path.join(aTask['wd'], rID + ".o.log"), 'a') as f:
                    f.write("\nCanceled with at least some jobs!")
                aTask['run'] = self.retryN + 1
            else:
                return True
        else:
            strF = os.path.join(aTask['wd'], rID + ".o.log")
            if os.path.isfile(strF):
                with open(strF, 'r') as f:
                    flines = f.readlines()
                if sum(bool(re.search('^DONE', line)) for line in flines) == aTask['task']:
                    aTask['run'] = -2
                    return False
                if aTask['task'] > 1:
                    fID = [int(line.split(": ")[-1]) for line in flines if line.startswith("DONE:")]
                    aID = [str(x) for x in (set(range(aTask['task'])) - set(fID))]
        if aTask['run'] > self.retryN:
            return True
        jID = self.run_slurm_job(aTask, rID, aID)
        if len(jID) < 3:
            print("Warning job submit error for %s" % os.path.join(aTask['wd'], rID + '.sh'))
            aTask['run'] += 1
        else:
            aTask['jID'] = jID
        return True
    
    def run_slurm_job(self, aTask, rID, aID=[], submit=True):
        if aTask['task'] == 1:
            subScript = self.slurmScript.replace("jID", os.path.basename(aTask['wd'])) \
                           .replace('jName', rID) \
                           .replace('wkPath', aTask['wd']) \
                           .replace('CoreNum', f"{aTask['cpu']*(1+aTask['run'])}") \
                           .replace('MEMFREE', f"{aTask['mem']*(1+aTask['run'])}") \
                           .replace('srcPath', self.strPipePath) \
                           .replace('RUMCMD', aTask['cmd'])
        else:
            subScript = self.slurmArray.replace("jID", os.path.basename(aTask['wd'])) \
                           .replace('jName', rID) \
                           .replace('wkPath', aTask['wd']) \
                           .replace('CoreNum', f"{aTask['cpu']*(1+aTask['run'])}") \
                           .replace('MEMFREE', f"{aTask['mem']*(1+aTask['run'])}") \
                           .replace('ARRAYN', f"{aTask['task']}") \
                           .replace('srcPath', self.strPipePath) \
                           .replace('RUMCMD', aTask['cmd'])
        fscript = os.path.join(aTask['wd'], rID + '.sh')
        with open(fscript, 'w') as f:
            f.write(subScript)
        if not submit:
            return ""
        # Wait until the script is written completely.
        while os.stat(fscript).st_size < 500:
            time.sleep(1)
        time.sleep(1)
        cmd = "sbatch " + ("" if len(aID) == 0 else "--array=%s " % ",".join(aID)) + f'{fscript}'
        print(f"submit: {os.path.basename(aTask['wd'])} {rID} by: {cmd}")
        result = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, text=True)
        return result.stdout.strip().split(" ")[-1]
    
    def get_slurm_job_status(self, jID):
        try:
            sacct_cmd = f'sacct --jobs={jID} --format=JobID,JobName%45,Partition,State,ExitCode,Start,End'
            result = subprocess.run(sacct_cmd, shell=True, stdout=subprocess.PIPE, text=True)
            df = pd.read_csv(StringIO(result.stdout), delim_whitespace=True)
        except Exception as e:
            print("Warnings with sacct: ", sacct_cmd)
            return 'RUNNING', []
        df = df[df.Partition == 'cpu']
        df['State'] = df.apply(lambda x: 'FAILED' if x['State'].startswith('OUT_OF_M') else x['State'], axis=1)
        State = df['State'].value_counts()
        if any(State.index.isin(['RUNNING','COMPLETING','CONFIGURING','RESIZING','PENDING'])):
            return 'RUNNING', []
        if any(State.index.isin(['CANCELLED','TIMEOUT','SUSPENDED','REVOKED','NODE_FAIL','PREEMPTED','RESIZING'])):
            return 'CANCELLED', []
        if all(State.index.isin(['COMPLETED'])):
            return 'COMPLETED', []
        if not all(State.index.str.startswith(('FAILED','COMPLETED'))):
            print(f"Warning unknown slurm states for {jID}: {','.join(State.index[~State.index.isin(['FAILED','COMPLETED'])].to_list())}")
            return 'RUNNING', []
        if df.shape[0] == 1:
            return 'FAILED', []
        df = df[df['State'] == 'FAILED']
        return 'FAILED', [x.split("_")[-1] for x in df['JobID'].to_list()]
    
    def resume_task(self):
        if self.specific_comparison is not None:
            return False
        comps = self.get_comparison()
        rJobs = []
        Resumed = False
        for one in comps:
            if not os.path.isfile(os.path.join(self.wdir, one, "run1.sh")):
                rJobs.append(one)
                continue
            strF = os.path.join(self.wdir, one, "run4.o.log")
            if os.path.isfile(strF):
                with open(strF, 'r') as f:
                    flines = f.readlines()
                if flines[-1].startswith("DONE"):
                    continue
            strF = os.path.join(self.wdir, one, "harm_run_4.e.log")
            if os.path.isfile(strF) and os.stat(strF).st_size == 0:
                continue
            rJobs.append(one)
            Resumed = True
        if len(rJobs) == 0:
            print("\n", self.wdir)
            print("*** All previous SplicingEvent jobs are completed in above folder!  ***")
            return True
        cmd = f"{self.strPipePath}/SpliceEvent_run {' '.join(sys.argv[1:])} -specific_comparison '{','.join(rJobs)}'"
        print(f"\nThe following {len(rJobs)} jobs will be resubmitted: ", "; ".join(rJobs))
        subprocess.run(f"nohup {cmd} &> {os.path.join(self.wdir, '.log_' + datetime.now().strftime('%Y%m%d_%H%M%S'))} &", shell=True)
        if Resumed:
            print("\n*** Resuming previous SplicingEvent jobs! Please check the status ***")
        else:
            print("\n*** Running SplicingEvent jobs! Please check the status ***")
        return True
    
    def process_comparisons(self, compList):
        cmd_dict = {}
        for comp in compList:
            cmd_dict[comp] = self.get_one_cmd(comp)
        if self.cmdonly:
            for k in cmd_dict:
                for k1 in sorted(cmd_dict[k].keys()):
                    self.run_slurm_job(cmd_dict[k][k1], k1, submit=False)
        else:
            while self.parallel_monitor(cmd_dict):
                time.sleep(60)
    
    def run_pipeline(self):
        if self.resume_task():
            return
        self.print_msg_init()
        compList = self.get_comparison()
        self.process_comparisons(compList)
        self.print_msg_power()


    @staticmethod
    def main(config_path):
        """Static main method that takes the config file path."""
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        pipeline = SpliceEventPipeline(config)
        pipeline.run_pipeline()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SpliceEvent pipeline using YAML config")
    parser.add_argument('-cfig', required=True, dest="cfig", help='YAML configuration file')
    args = parser.parse_args()
    start_time = time.time()
    SpliceEventPipeline.main(args.cfig)
    print("\n--- total time passed %s seconds ---" % (time.time() - start_time))