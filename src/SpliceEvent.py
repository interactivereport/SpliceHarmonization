import argparse,os,shutil,glob,ast,time,configparser,yaml,pickle,re,time,sys,subprocess
import pandas as pd
from datetime import datetime
from io import StringIO
from tabulate import tabulate
from slurmUtils import slurmScript, slurmArray

class SpliceEventPipeline:
    def __init__(self):
        self.strPipePath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.GTFpath = '/edgehpc/dept/compbio/reference/DNAnexus_references'
        self.retryN = 2
        self.g_config = None
    
    def exit_with_msg(self, msg=""):
        if len(msg)>3:
            print(msg)
        sys.exit(1)

    def getConfig(self, strConfig):
        with open(strConfig,'r') as f:
            self.g_config = yaml.safe_load(f)
    
    def msg_init(self):
        print("\n\n*****",datetime.now().strftime("%Y-%m-%d %H:%M:%S"),"*****")
        if os.path.isdir(os.path.join(self.strPipePath,".git")):
            gitConfig = configparser.ConfigParser()
            tmp = gitConfig.read(os.path.join(self.strPipePath,".git","config"))
            url = gitConfig['remote "origin"']['url']
            gitLog = pd.read_csv(os.path.join(self.strPipePath,".git","logs","HEAD"),sep="\t",header=None)
            gitLog = gitLog.iloc[-1,0].split(" ")
            rInfo=[]
            rInfo.append(f"###########\n## SpliceEvent: {url}\n")
            rInfo.append(f"## Run Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            rInfo.append(f"## Pipeline Path: {self.strPipePath}\n")
            rInfo.append(f"## Pipeline Date: {datetime.fromtimestamp(int(gitLog[-2])).strftime('%Y-%m-%d %H:%M:%S')} {gitLog[-1]}\n")
            rInfo.append(f"## git HEAD: {gitLog[1]}\n###########\n\n")
            rInfo.append(f"{self.strPipePath}/SpliceEvent_run {' '.join(sys.argv[3:])}\n")
            print(rInfo)
            print(''.join(rInfo))
            with open(os.path.join(self.g_config['wdir'],'_cmdline_'+datetime.now().strftime('%Y%m%d.%H%M%S')),'w') as f:
                f.writelines(rInfo)
            #print("## Run Date: ",datetime.now())
            #print("## Pipeline Path: %s"%self.strPipePath)
            #print("## Pipeline Date: %s %s"%(datetime.fromtimestamp(int(gitLog[-2])).strftime('%Y-%m-%d %H:%M:%S'),gitLog[-1]))
            #print("## git HEAD: %s\n###########\n"%gitLog[1])

    def get_comparison(self):
        stringtie_path = os.path.join(self.g_config['indir'], "stringtie")
        if not os.path.isdir(stringtie_path):
            self.exit_with_msg("Missing stringtie folder: %s"%stringtie_path)
        comp = [re.sub('.combined.gtf$','', os.path.basename(_)) for _ in glob.glob(os.path.join(stringtie_path,"*.combined.gtf"))]
        if self.g_config['specific_comparison'] is not None:
            speComp = self.g_config['specific_comparison'].split(',')
            comp = [_ for _ in comp if any(s in _ for s in speComp)]
        if len(comp) == 0:
            self.exit_with_msg("No specific comparison found from stringtie")
        return(comp)

    def process_comparison(self,cList):
        cmd = {}
        for one in cList:
            cmd[one] = self.get_one_cmd(one)
        if self.g_config.cmdonly:
            for k in cmd:
                for k1 in sorted(cmd[k].keys()):
                    self.run_slurm_job(cmd[k][k1],k1,submit=False)
        else:
            while self.parallel_monitor(cmd,self.g_config['wdir']):
                time.sleep(60)

    def get_one_cmd(self,oneComp, args):
        data_folder = self.g_config['indir']
        bam_folder = self.g_config['bam_path']
        wd_folder = self.g_config['wdir']
        reference = self.g_config['genome_gtf']
        specific_comparison_str = self.g_config['specific_comparison']
        sample_sheet_file = self.g_config['sample_info']
        file_prefix = os.path.join(wd_folder,oneComp)
        os.makedirs(os.path.join(file_prefix,'junction_prep'), exist_ok=True)

        comparisonRef_name = oneComp.split('_vs_')[-1]
        rmats_input = os.path.join(data_folder,'rmats',oneComp+'_*JCEC.txt' )
        leafcutter_input = os.path.join(data_folder,'leafcutter',oneComp+'*_cluster_significance.txt' )
        majiq_input = os.path.join(data_folder,'majiq',oneComp+'*voila.tsv' )
        rmats_output= os.path.join(file_prefix,'junction_prep','rmats_junction_prep.csv')
        leafcutter_output= os.path.join(file_prefix,'junction_prep','leafcutter_junction_prep.csv')
        majiq_output= os.path.join(file_prefix,'junction_prep','majiq_junction_prep.csv')
        
        junction_input= os.path.join(file_prefix,'junction_prep','{method}_junction_prep.csv')
        sj_output = os.path.join(file_prefix,'junction_prep','sj.csv')
        
        annotation= os.path.join(data_folder,'stringtie',oneComp+'.combined.gtf')
        harm_prep_dir= os.path.join(file_prefix,'harm')
        mxe_prep_dir = os.path.join(file_prefix,'harm','mxe')
        harm_post_input =  os.path.join(harm_prep_dir,'harm*.csv')
        mxe_post_input = os.path.join(mxe_prep_dir,'harm*.csv')
        harm_post_dir= os.path.join(file_prefix,'harm','out')
        mxe_post_dir= os.path.join(file_prefix,'harm','mxe','out')
        relabel = os.path.join(harm_prep_dir,'exon_relabel.csv')
        majiq_cutoff = os.path.join(harm_prep_dir,'majiq_cutoff_value.csv')

        final_output= os.path.join(file_prefix,'out')
        final_mxe_output= os.path.join(file_prefix,'out','mxe')
        final_allgene_output=  os.path.join(file_prefix,'out','all_gene')
        all_gene_input =  os.path.join(final_allgene_output,'*.csv')
        all_junction_output =  os.path.join(file_prefix,'out','all_gene_junction')
        
        os.makedirs(harm_post_dir, exist_ok=True)
        os.makedirs(mxe_post_dir, exist_ok=True)
        os.makedirs(final_mxe_output, exist_ok=True)
        os.makedirs(final_allgene_output, exist_ok=True)
        os.makedirs(all_junction_output, exist_ok=True)
        
        bam_pattern = "Aligned.sortedByCoord.out.bam"
        if args.junction_filter and sample_sheet_file:
            sample_sheet = pd.read_csv(sample_sheet_file,sep="\t")
            sample_sheet["SampleName"]=sample_sheet["SampleName"].astype(str)
            sample_sheet["Groups"]=sample_sheet["Groups"].astype(str)
            #Iterate through rows of grouping file to create individual group files per analysis
            alt_rows = sample_sheet[sample_sheet['Groups'] == file_prefix.split('_vs_')[0]].index
            print(file_prefix.split('_vs_')[0])
            alt_bams = sample_sheet.SampleName[alt_rows] + '.' + bam_pattern
            print(alt_bams)
            ref_rows = sample_sheet[sample_sheet['Groups'] == file_prefix.split('_vs_')[1]].index
            ref_bams = sample_sheet.SampleName[ref_rows] + '.' + bam_pattern
        if args.junction_filter and bam_folder:
            altbam_files = [os.path.join(bam_folder,x) for x in alt_bams.to_list()]
            refbam_files = [os.path.join(bam_folder,x) for x in ref_bams.to_list()]
        inputjson=os.path.join(final_output,'libdepth_counts.json')
        
        allCMD = {}
        
        # run 1
        strCMD = []
        strCMD.append(f"Rscript %s -i '{rmats_input}' -o {rmats_output} -n {comparisonRef_name}"%os.path.join(self.strPipePath,"src/SpliceEvent","rmats_junction_prep.R"))
        strCMD.append(f"Rscript %s -i '{leafcutter_input}' -o {leafcutter_output}"%os.path.join(self.strPipePath,"src/SpliceEvent","leafcutter_junction_prep.R"))
        strCMD.append(f"Rscript %s -i '{majiq_input}' -o {majiq_output}"%os.path.join(self.strPipePath,"src/SpliceEvent","majiq_junction_prep.R"))
        strCMD.append(f"Rscript %s -a {annotation} -r {reference} -i {junction_input} -o {sj_output}"%os.path.join(self.strPipePath,"src/SpliceEvent","annotation_SJ.R"))
        strCMD.append(f"python -u %s -a {annotation} -r {reference} -sj {sj_output} -o {harm_prep_dir}"%os.path.join(self.strPipePath,"src/SpliceEvent","harm_prep_remove_mapping_add_exons.py") + (" -addexons %s"%args.addexons if args.addexons else ""))
        allCMD['run1'] = {'cmd':'\n'.join(strCMD),'wd':file_prefix,'mem':96,'cpu':8,'task':1,'run':0}
        
        # run2
        strCMD.clear()
        strCMD.append(f"combined_files=($(ls {harm_post_input}) $(ls {mxe_post_input}))")
        strCMD.append("input_csv=${combined_files[$SLURM_ARRAY_TASK_ID-1]}")
        strCMD.append(f"python -u %s -exonlabel {relabel} -infile $input_csv -outdir {harm_post_dir} -majiqcutofffile {majiq_cutoff}"%os.path.join(self.strPipePath,"src/SpliceEvent","harm_assign_job_array.py"))
        allCMD['run2'] = {'cmd':'\n'.join(strCMD),'wd':file_prefix,'mem':16,'cpu':1,'task':80,'run':0}

        # run 3
        strCMD.clear()
        strCMD.append(f"python -u %s -indir {harm_post_dir} -outdir {final_output} -allgene {final_allgene_output} -majiq_score_cutoff {args.majiq_cutoff_val}"%os.path.join(self.strPipePath,"src/SpliceEvent","post_process.py"))
        strCMD.append(f"python -u %s -indir {mxe_post_dir} -outdir {final_mxe_output}"%os.path.join(self.strPipePath,"src/SpliceEvent","post_process_4mxe.py"))
        if args.junction_filter:
            strCMD.append(f"python -u %s -outdir {final_output} -Altbam {altbam_files} -Refbam {refbam_files}"%os.path.join(self.strPipePath,"src/SpliceEvent","libdepth_count.py"))
        allCMD['run3'] = {'cmd':'\n'.join(strCMD),'wd':file_prefix,'mem':96,'cpu':8,'task':1,'run':0}
        
        # run 4
        strCMD.clear()
        if args.junction_filter:
            strCMD.append(f"combined_files=($(ls {all_gene_input}) %s)"%os.path.join(final_mxe_output,"all_cleaned.csv"))
            strCMD.append("files_per_task=$(( (${#combined_files[@]} + SLURM_ARRAY_TASK_MAX - 1) / SLURM_ARRAY_TASK_MAX ))")
            strCMD.append("start_idx=$(( (SLURM_ARRAY_TASK_ID - 1) * files_per_task ))")
            strCMD.append("end_idx=$(( start_idx + files_per_task - 1 ))")
            strCMD.append("if [ $end_idx -ge ${#combined_files[@]} ]; then\n\tend_idx=$((${#combined_files[@]} - 1))\nfi")
            strCMD.append('for i in $(seq $start_idx $end_idx); do\n\tinput_csv="${combined_files[$i]}"')
            strCMD.append(f'\tpython -u %s -infile $input_csv -outdir {all_junction_output} -Altbam "{altbam_files}" -Refbam "{refbam_files}" -inputjson {inputjson} -FCthreshold {args.junctionFC} -maxcount {args.junctionMAX}\ndone')
            allCMD['run4'] = {'cmd':'\n'.join(strCMD),'wd':file_prefix,'mem':8,'cpu':1,'task':50,'run':0}
            strCMD.clear()
            strCMD.append(f"python -u %s -indir {all_junction_output} -mxeindir {final_mxe_output} -outdir {final_mxe_output}"%os.path.join(self.strPipePath,"src/SpliceEvent","post_process_MXE_junction.py"))
            allCMD['run5'] = {'cmd':'\n'.join(strCMD),'wd':file_prefix,'mem':16,'cpu':1,'task':1,'run':0}
        else:
            strCMD.append(f"python -u %s -indir {final_allgene_output} -mxeindir {final_mxe_output} -outdir {final_mxe_output}"%os.path.join(self.strPipePath,"src/SpliceEvent","post_process_MXE.py"))
            allCMD['run4'] = {'cmd':'\n'.join(strCMD),'wd':file_prefix,'mem':96,'cpu':8,'task':1,'run':0}
            
        
        return allCMD

    def parallel_monitor(self, cmdList,wd):
        stat = {}
        for k in cmdList:
            for k1 in sorted(cmdList[k].keys()):
                if cmdList[k][k1]['run']<0:
                    continue
                if self.check_task(cmdList[k][k1],k1):
                    break
        #status
        stat = {}
        running = False
        for k in cmdList:
            wait = False
            for k1 in sorted(cmdList[k].keys()):
                if not k1 in stat:
                    stat[k1] = {}
                if wait:
                    stat[k1][k] = "Waiting"
                    continue
                match cmdList[k][k1]['run']:
                    case x if x==-2:
                        stat[k1][k] = "Resumed"
                    case x if x==-1:
                        stat[k1][k] = "Completed"
                    case x if x>self.retryN:
                        stat[k1][k] = "Failed"
                    case x:
                        stat[k1][k] = f"Running_{cmdList[k][k1]['jID'] if 'jID' in cmdList[k][k1] else ''}_{x}"
                        running=True
                if stat[k1][k].startswith('Running') or stat[k1][k].startswith('Failed'):
                    wait=True
        with open(os.path.join(wd,'status.log'),'w') as f:
            f.write(tabulate(pd.DataFrame(stat), headers='keys',tablefmt='plain')+"\n\n")
        #pd.DataFrame(stat).to_csv(os.path.join(wd,'status.csv'))
        return running

    def check_task(self, aTask,rID):
        aID = []
        if 'jID' in aTask.keys():
            jStat,aID = self.get_slurm_job_status(aTask['jID'])
            if jStat=="COMPLETED":
                aTask['run'] = -1
                return False
            elif jStat=="RUNNING":
                return True
            elif jStat=="FAILED":
                if os.path.isfile(os.path.join(aTask['wd'],rID+".e.log")):
                    shutil.move(os.path.join(aTask['wd'],rID+".e.log"),
                        os.path.join(aTask['wd'],rID+"_%d.e.log"%aTask['run']))
                with open(os.path.join(aTask['wd'],rID+".o.log"),'a') as f:
                    f.write("\nFailed with at least some jobs!")
                aTask['run'] += 1
            elif jStat=='CANCELLED':
                with open(os.path.join(aTask['wd'],rID+".o.log"),'a') as f:
                    f.write("\nCanceled with at least some jobs!")
                aTask['run'] = self.retryN+1
            else:#other status
                return True
        else:
            strF = os.path.join(aTask['wd'],rID+".o.log")
            if os.path.isfile(strF):
                with open(strF,'r') as f:
                    flines = f.readlines()
                if sum(bool(re.search('^DONE',_)) for _ in flines)==aTask['task']:
                    aTask['run'] = -2
                    return False
                if aTask['task']>1:
                    fID = [int(_.split(": ")[-1]) for _ in flines if _.startswith("DONE:")]
                    aID = [str(_) for _ in (set(range(aTask['task']))-set(fID))]
        if aTask['run']>self.retryN:
            return True
        jID = self.run_slurm_job(aTask,rID,aID)
        if len(jID)<3:
            print("Warning job submit error for %s"%os.path.join(aTask['wd'],rID+'.sh'))
            aTask['run'] += 1
        else:
            aTask['jID'] = jID
        return True

    def run_slurm_job(self, aTask,rID,aID=[],submit=True):
        if aTask['task']==1:
            subScript = (slurmScript.
                replace("jID",os.path.basename(aTask['wd'])).
                replace('jName',rID).
                replace('wkPath',aTask['wd']).
                replace('CoreNum',f"{aTask['cpu']*(1+aTask['run'])}").
                replace('MEMFREE',f"{aTask['mem']*(1+aTask['run'])}").
                replace('srcPath',self.strPipePath).
                replace('RUMCMD',aTask['cmd']))
        else:
            subScript = (slurmArray.
                replace("jID",os.path.basename(aTask['wd'])).
                replace('jName',rID).
                replace('wkPath',aTask['wd']).
                replace('CoreNum',f"{aTask['cpu']*(1+aTask['run'])}").
                replace('MEMFREE',f"{aTask['mem']*(1+aTask['run'])}").
                replace('ARRAYN',f"{aTask['task']}").
                replace('srcPath',self.strPipePath).
                replace('RUMCMD',aTask['cmd']))
        fscript = os.path.join(aTask['wd'],rID+'.sh')
        with open(fscript,'w') as f:
            f.write(subScript)
        if not submit:
            return
        # wait for file writing completed before submit to cluster
        while os.stat(fscript).st_size<500:
            time.sleep(1)
        time.sleep(1)
        cmd = "sbatch "+("" if len(aID)==0 else "--array=%s "%",".join(aID))+f'{fscript}'
        print(f"submit: {os.path.basename(aTask['wd'])} {rID} by: {cmd}")
        return subprocess.run(cmd,shell=True,stdout=subprocess.PIPE, text=True).stdout.strip().split(" ")[-1]
        
    def get_slurm_job_status(self, jID):
        try:
            df = pd.read_csv(StringIO(subprocess.run(f'sacct --jobs={jID} --format=JobID,JobName%45,Partition,State,ExitCode,Start,End',
                shell=True,stdout=subprocess.PIPE, text=True).stdout),delim_whitespace=True)
        except:
            print("Warnings with sacct: ",f'sacct --jobs={jID} --format=JobID,JobName%45,Partition,State,ExitCode,Start,End')
            return 'RUNNING',[]
        df = df[df.Partition=='cpu']
        df.State = df.apply(lambda x: 'FAILED' if x['State'].startswith(('OUT_OF_M')) else x['State'],axis=1)
        State = df.State.value_counts()
        if any(State.index.isin(['RUNNING','COMPLETING','CONFIGURING','RESIZING','PENDING'])):
            return 'RUNNING',[]
        if any(State.index.isin(['CANCELLED','TIMEOUT','SUSPENDED','REVOKED','NODE_FAIL','PREEMPTED','RESIZING'])):
            return 'CANCELLED',[]
        if all(State.index.isin(['COMPLETED'])):
            return 'COMPLETED',[]
        if not all(State.index.str.startswith(('FAILED','COMPLETED'))):
            print(f"Warning unknown slurm states for {jID}: {','.join(State.index[~State.index.isin(['FAILED','COMPLETED'])].to_list())}")
            return 'RUNNING',[]
        if df.shape[0]==1:
            return 'FAILED',[]
        df = df[df.State=='FAILED']
        return 'FAILED',[_.split("_")[-1] for _ in df.JobID.to_list()]

    def resume_task(self):
        #print(args)
        if self.specific_comparison is not None:
            return False
        comp = self.get_comparison()
        rJobs = []
        Resumed=False
        # print(self.g_args.wdir)
        for one in comp:
            if not os.path.isfile(os.path.join(self.g_config['wdir'],one,"run1.sh")):
                rJobs.append(one)
                continue
            strF = os.path.join(self.g_config['wdir'],one,"run4.o.log")
            if os.path.isfile(strF):
                with open(strF,'r') as f:
                    flines = f.readlines()
                if flines[-1].startswith("DONE"):
                    continue
            # old pipeline results
            strF = os.path.join(self.g_config['wdir'],one,"harm_run_4.e.log")
            if os.path.isfile(strF) and os.stat(strF).st_size==0:
                continue
            rJobs.append(one)
            Resumed=True
        if len(rJobs)==0:
            print("\n",self.g_config['wdir'])
            print("*** All previous SplicingEvent jobs are completed in above folder!  ***")
            return True
        cmd = (f"{self.strPipePath}/SpliceEvent.py " + " ".join(sys.argv[1]) +" -specific_comparison '%s'" % (",".join(rJobs)))
        print(cmd)
        print("\nThe following %d jobs will be resubmited: "%len(rJobs),"; ".join(rJobs))
        
        a = subprocess.run(f"nohup {cmd} &> {os.path.join(self.g_config['wdir'],'.log_'+datetime.now().strftime('%Y%m%d_%H%M%S'))} &",shell=True)
        if Resumed:
            print("\n*** Resuming prevous SplicingEvent jobs! Please check the status ***")
        else:
            print("\n*** Runing SplicingEvent jobs! Please check the status ***")
        return True


    # def run(self):
    #     if not os.path.isdir(args.indir):
    #         self.exit_with_msg('Missing fastr splicing folder: %s'%args.indir)
    #     if not os.path.isdir(os.path.dirname(args.wdir)):
    #         self.exit_with_msg('Missing output folder: %s'%os.path.dirname(args.wdir))
    #     os.makedirs(args.wdir,exist_ok=True)
    #     if args.ref is None:
    #         strAna = glob.glob(os.path.join(args.indir,'analysis-*'))[0]
    #         if os.path.isfile(strAna):
    #             with open(strAna,'r') as f:
    #                 aInfo = yaml.safe_load(f)
    #             version = aInfo["FASTR_REFERENCE_VERSION"]
    #             args.ref = os.path.join(self.GTFpath,'rnaseq',
    #                 aInfo["FASTR_REFERENCE_SPECIES"],
    #                 version,"%s.transcript.gtf.gz"%version)
    #     if args.ref is None or not os.path.isfile(args.ref):
    #         self.exit_with_msg('Missing reference GTF: %s'%args.ref)
    #     if args.junction_filter:
    #         if args.samplesheet is None or args.bamindir is None or not os.path.isfile(args.samplesheet) or not os.path.isdir(args.bamindir):
    #             self.exit_with_msg(f'junction_filter (True) requires samplesheet ({args.samplesheet}) and bamindir ({args.bamindir})')

    #     if self.resume_task():
    #         return
    
        # self.msg_init(args.wdir)
        # compList = self.get_comparison()
        # self.process_comparison(self.args,compList)


def main(strConfig):
    # Create an instance of the class
   
    # self.getConfig(strConfig)
    pipeline = SpliceEventPipeline() 
    pipeline.getConfig(strConfig) 
    if pipeline.resume_task():
        return
    
    pipeline.msg_init()
    compList = pipeline.get_comparison()
    pipeline.process_comparison(compList)
# Standalone fallback: If this module is run directly, parse global args.
if __name__ == '__main__':

    main()
