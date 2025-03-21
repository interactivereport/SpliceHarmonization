import os,sys,yaml,subprocess,re,glob,shutil,time
import pandas as pd
from datetime import datetime

class SplicePrepPipeline:
    def __init__(self, config):
        self.g_pipePath = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.g_debug = False
        self.g_config = config
        self.g_sysConfig = None
        self.g_script=None

    def msgError(self, msg):
        print('\n'.join(msg))
        sys.exit(1)

    def getSysConfig(self):
        with open(os.path.join(self.g_pipePath,'src/SplicePrep','sys.yml'),'r') as f:
            self.g_sysConfig = yaml.safe_load(f)

    @staticmethod
    def getConfig(strConfig):
        with open(strConfig,'r') as f:
            config = yaml.safe_load(f)
        return config

    def checkConfig(self):
        errorMsg = []
        if not os.path.isfile(self.g_config['sample_info']):
            errorMsg.append('Missing the sample info file: %s'%self.g_config['sample_info'])
        if not os.path.isfile(self.g_config['compare_info']):
            errorMsg.append('Missing the compare info file: %s'%self.g_config['compare_info'])
        if not os.path.isfile(self.g_config['genome_fa']):
            errorMsg.append('Missing the genome fasta file: %s'%self.g_config['genome_fa'])
        if not os.path.isfile(self.g_config['genome_gtf']):
            errorMsg.append('Missing the genome GTF file: %s'%self.g_config['genome_gtf'])
        if not os.path.isdir(self.g_config['bam_path']):
            errorMsg.append('Missing the BAM file path: %s'%self.g_config['bam_path'])
        if self.g_config['genome_fa'].endswith('gz'):
            errorMsg.append('genome fasta file is required to be unzipped: %s'%self.g_config['genome_fa'])
        if self.g_config['genome_gtf'].endswith('gz'):
            errorMsg.append('genome GTF file is required to be unzipped: %s'%self.g_config['genome_gtf'])
        if len(errorMsg)>0:
            self.msgError(errorMsg)
        self.getSampleGroup()
        if not 'Sample_Name' in self.g_sInfo.columns:
            errorMsg.append('Missing Sample_Name column in the sample file: %s'%self.g_config['sample_info'])
        if not set(['Group', 'Alt', 'Ref']).issubset(self.g_cInfo.columns):
            errorMsg.append('Missing some required columns (Group, Alt, Ref) in the comparison file: %s'%self.g_config['compare_info'])
        if len(errorMsg)>0:
            self.msgError(errorMsg)
        if not set(self.g_cInfo.Group).issubset(self.g_sInfo.columns):
            errorMsg.append('Missing some group columns from the sample file')
        for one in self.g_sInfo.Sample_Name:
            strBAM = os.path.join(self.g_config['bam_path'],one+"."+self.g_config['bam_suffix'])
            if not os.path.isfile(strBAM):
                errorMsg.append('Missing sample BAM files for %s: %s'%(one,strBAM))
        for one in self.g_cInfo.Group.unique():
            if not set(self.g_cInfo[self.g_cInfo.Group==one][['Alt','Ref']].to_numpy().flatten()).issubset(self.g_sInfo[one]):
                errorMsg.append('Missing some values in Group (%s) from the sample file'%one)
        if len(errorMsg)>0:
            self.msgError(errorMsg)

    def getSampleGroup(self):
        self.g_sInfo = pd.read_csv(self.g_config['sample_info'], sep='\t')
        self.g_sInfo['Sample_Name'] = self.g_sInfo['Sample_Name'].astype(str)
        self.g_cInfo = pd.read_csv(self.g_config['compare_info'], sep='\t')

    def getSbatch(self):
        with open(os.path.join(self.g_pipePath,'src/SplicePrep','subjob.sh'),'r') as f:
            self.g_script = f.read()

    def getOneSbatch(self, jName,wkPath,cmd,CoreNum=4,MEMFREE=None):
        if MEMFREE is None:
            MEMFREE = 16*CoreNum
        
        oneSub = self.g_script.format(jName=jName,
            wkPath=wkPath,
            CoreNum=CoreNum,
            MEMFREE=MEMFREE,
            strCMD='\n'.join(cmd))
        oneF = os.path.join(wkPath,jName+'.sh')
        with open(oneF,'w') as f:
            f.write(oneSub)
        return oneF

    # def getOneSbatch(self, jName, output_path, cmd, n_tasks=4, mem_free=None, time_limit="72:0:0"):
    #     script = create_sbatch_script(jName, output_path, cmd, n_tasks, mem_free, time_limit)
    #     script_path = write_sbatch_script(script, output_path, jName)
    #     return script_path
    
    def saveGrpFile(self, one,strF):
        with open(strF,'w') as f:
            for _ in self.g_sInfo[self.g_sInfo[one.Group]==one.Ref].Sample_Name:
                f.write("%s\t%s\n"%(_,one.Ref))#,self.g_config['bam_suffix']
            for _ in self.g_sInfo[self.g_sInfo[one.Group]==one.Alt].Sample_Name:
                f.write("%s\t%s\n"%(_,one.Alt))#,self.g_config['bam_suffix']

    def initRun(self, job):
        if not job in self.g_config['run']:
            print("Skip %s!"%job)
            return None,None,None
        print("Submitting %s jobs ..."%job)
        strOut = os.path.join(self.g_config['output_path'],'%s_output'%job)
        strTmp = os.path.join(strOut,'tmp')
        strSbatch = os.path.join(strOut,'sbatch')
        os.makedirs (strTmp,exist_ok=True)
        os.makedirs (strSbatch,exist_ok=True)
        return strOut,strTmp,strSbatch
    
    def getBAMindex(self):
        print("Submitting the bam index")
        strTmp = os.path.join(self.g_config['output_path'],'BAM_index_sbatch')
        os.makedirs(strTmp,exist_ok=True)
        jobIDs = []
        for sName in self.g_sInfo.Sample_Name:
            strBAM = os.path.join(self.g_config['bam_path'],sName+"."+self.g_config['bam_suffix'])
            if not os.path.isfile(strBAM+'.bai'):
                cmd = ['module load SAMtools']
                cmd.append('samtools index -@ %d %s'%(4,strBAM))
                if self.g_debug:
                    print('\tix_%s'%sName)
                else:
                    jobIDs.append(subprocess.run("sbatch %s"%self.getOneSbatch('ix_%s'%sName,strTmp,cmd,4),
                        shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        return jobIDs

    def stringtie_run(self):
        strOut,strTmp,strSbatch = self.initRun('stringtie')
        if strOut is None:
            return None
        jIDs = {}
        for i in range(self.g_cInfo.shape[0]):
            for g in ['Alt','Ref']:
                strJob = self.g_cInfo.Group[i]+"_"+self.g_cInfo[g][i]
                oneFinal = os.path.join(strTmp,strJob+'.gtf')
                if os.path.isfile(oneFinal) or strJob in jIDs.keys():
                    continue
                oneBatch = self.getOneSbatch('ST_'+strJob,strSbatch,
                    self.stringtie_one(self.g_cInfo.Group[i],self.g_cInfo[g][i],strTmp,oneFinal))
                if self.g_debug:
                    print("\t%s"%strJob)
                else:
                    jIDs[strJob] = subprocess.run("sbatch %s"%oneBatch,shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip()
        allJobs=[]
        job_dependency = '--dependency=afterok:%s'%','.join(jIDs.values()) if len(jIDs)>0 else ""
        for i in range(self.g_cInfo.shape[0]):
            strCom = self.g_cInfo.Group[i]+'_g_'+ self.g_cInfo.Alt[i]+"_vs_"+ self.g_cInfo.Ref[i]
            oneTmp = os.path.join(strTmp,strCom)
            os.makedirs(oneTmp,exist_ok=True)
            strGFF_ref = os.path.join(oneTmp, self.g_cInfo.Group[i]+"_"+ self.g_cInfo.Ref[i]+'.gtf') #check! might not needed to copy if gffcompare doesn't change the file
            strGFF_alt = os.path.join(oneTmp, self.g_cInfo.Group[i]+"_"+ self.g_cInfo.Alt[i]+'.gtf')
            prefix = os.path.join(oneTmp,self.g_cInfo.Alt[i]+"_vs_"+self.g_cInfo.Ref[i])
            if os.path.isfile(prefix+'.combined.gtf'):
                continue
            cmd = ['module purge &> /dev/null','module load GffCompare']
            cmd.append("cp %s %s"%(os.path.join(strTmp, self.g_cInfo.Group[i]+"_"+self.g_cInfo.Ref[i]+'.gtf'),strGFF_ref))
            cmd.append("cp %s %s"%(os.path.join(strTmp, self.g_cInfo.Group[i]+"_"+self.g_cInfo.Alt[i]+'.gtf'),strGFF_alt))
            cmd.append("gffcompare -o %s -r %s -R -V %s %s"%(prefix,self.g_config['genome_gtf'],strGFF_ref,strGFF_alt))
            oneBatch = self.getOneSbatch('ST_'+strCom,strSbatch,cmd,CoreNum=2)
            if self.g_debug:
                print("\t%s"%strCom)
            else:
                allJobs.append(subprocess.run("sbatch %s %s"%(job_dependency,oneBatch),#--export=NONE
                    shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        return allJobs
    
    def stringtie_one(self, k,v,strTmp,strFinal):
        strBAM = [os.path.join(self.g_config['bam_path'],_+"."+self.g_config['bam_suffix']) for _ in self.g_sInfo[self.g_sInfo[k]==v]['Sample_Name']]
        strG = os.path.join(strTmp,"%s_%s"%(k,v))
        library_type = "" if self.g_config['strand']==0 else ("--rf" if self.g_config['strand']==1 else "--fr")
        cmd = ["module load StringTie",'module load SAMtools']
        cmd.append('samtools merge -f -@ 4 %s.bam %s'%(strG,' '.join(strBAM)))
        cmd.append('samtools view -h %s.bam | gawk -v strType=2 -f %s | samtools view -bS > %s.xs.bam'%(strG,os.path.join(self.g_pipePath,'src/SplicePrep','tagXSstrandedData.awk'),strG))
        cmd.append("stringtie %s.xs.bam -o %s -p 16 %s -c %d -j %d -f %f"%(
            strG,strFinal,library_type,self.g_config['stringtie_min_cov'],self.g_config['stringtie_min_jxn_cov'],self.g_config['stringtie_min_isoform_frac']))
        return cmd

    def leafcutter_run(self, index_jobID):
        strOut,strTmp,strSbatch = self.initRun('leafcutter')
        if strOut is None:
            return None
        strGTF2leafcutter = os.path.join(strTmp,re.sub(".gtf$","",os.path.basename(self.g_config['genome_gtf'])))#+"_all_exons.txt.gz"
        # gtf to leafcutter
        jobIDs =  []
        if not os.path.isfile(strGTF2leafcutter+"_all_exons.txt.gz"):
            cmd = []
            cmd.append('%s --output %s %s'%(os.path.join(self.g_pipePath,'src/SplicePrep','leafcutter_gtf2leafcutter.pl'),strGTF2leafcutter,self.g_config['genome_gtf']))
            oneBatch = self.getOneSbatch('LC_GTF',strSbatch,cmd,CoreNum=2)
            if self.g_debug:
                print("\tgenome GTF to leafcutter")
            else:
                jobIDs.append(subprocess.run("sbatch %s"%oneBatch,shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        # bam to junc
        index_dependency = '--dependency=afterok:%s'%','.join(index_jobID) if len(index_jobID)>0 else ''
        for i in range(self.g_sInfo.shape[0]):
            oneScript = self.leafcutter_bam2junc(self.g_sInfo.Sample_Name[i],strTmp,strSbatch)
            if oneScript is None:
                continue
            else:
                if self.g_debug:
                    print("\tjunc:",self.g_sInfo.Sample_Name[i])
                else:
                    jobIDs.append(subprocess.run("sbatch %s %s"%(index_dependency,oneScript),
                    shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        #juct to clust
        job_dependency = '--dependency=afterok:%s'%','.join(jobIDs) if len(jobIDs)>0 else ''
        allJobs=[]
        for row in self.g_cInfo.itertuples():
            oneScript = self.leafcutter_diffSplicing(row,strTmp,strGTF2leafcutter,strSbatch)
            if oneScript is None:
                continue
            if self.g_debug:
                print("\t%s"%os.path.basename(oneScript))
            else:
                allJobs.append(subprocess.run("sbatch %s %s"%(job_dependency,oneScript),#--export=NONE 
                    shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        return allJobs
    def leafcutter_bam2junc(self, sName,strTmp,strSbatch):
        #print("\tBam to Junc: ",sName)
        strBAM = os.path.join(self.g_config['bam_path'],sName+"."+self.g_config['bam_suffix'])
        strJunc = os.path.join(strTmp,sName+".junc")
        if os.path.isfile(strJunc):
            return None
        cmd = ['module load RegTools']
        cmd.append('regtools junctions extract -a %d -m %d -M %d -s %d %s -o %s'%(
            self.g_config['leafcutter_anchor_len'],self.g_config['leafcutter_min_intron_len'],self.g_config['leafcutter_max_intron_len'],
            self.g_config['strand'],strBAM,strJunc
            ))
        return self.getOneSbatch('LC_'+sName,strSbatch,cmd,4)
    
    def leafcutter_diffSplicing(self, one,strTmp,strGTF2leafcutter,strSbatch):
        strCom = one.Group+"_g_"+one.Alt+"_vs_"+one.Ref
        #print("\t",strCom)
        strJunc = [os.path.join(strTmp,_+".junc") for _ in self.g_sInfo[self.g_sInfo[one.Group].isin(one)].Sample_Name]
        oneTmp = os.path.join(strTmp,strCom)
        os.makedirs(oneTmp,exist_ok=True)
        grpF = os.path.join(oneTmp,'grp_files.txt')
        listF = os.path.join(oneTmp,'junc_files.txt')
        junc2clust_prefix = os.path.join(oneTmp,'junct2clust')
        leafcutterDS_prefix = os.path.join(oneTmp,one.Alt+"_vs_"+one.Ref+'_ds_res')
        leafcutterDS_plot = os.path.join(oneTmp,'leafcutter.ds.pdf')
        leafcutterDS_rdata = os.path.join(oneTmp,one.Alt+"_vs_"+one.Ref+'_leafviz.Rdata')
        if os.path.isfile(leafcutterDS_rdata):
            return None
        allF=[]
        cmd = ['shopt -s expand_aliases','module load Python/2.7.18-GCCcore-12.2.0-bare']
        for _ in strJunc:
            allF.append(os.path.join(oneTmp,os.path.basename(_)))
            if self.g_config['strand']==0:
                cmd.append("cp %s %s"%(_,allF[-1]))
            else:
                cmd.append("grep -v '?' %s > %s"%(_,allF[-1]))
        #allF.append("")
        with open(listF,'w') as f:
            a = f.write('\n'.join(allF))
        #cmd.append("cd %s"%os.path.dirname(listF))
        cmd.append("python %s --juncfiles %s --minclureads %d --mincluratio %f --outprefix %s --maxintronlen %d --rundir %s"%(
            os.path.join(self.g_pipePath,'src/SplicePrep','leafcutter_cluster_regtools.py'),
            listF,
            self.g_config['leafcutter_min_cluster_reads'],
            self.g_config['leafcutter_min_cluster_ratio'],
            os.path.basename(junc2clust_prefix),
            self.g_config['leafcutter_max_intron_len'],
            oneTmp
        ))
        cmd.append("module load Leafcutter")
        self.saveGrpFile(one,grpF)
        cmd.append("Rscript %s --num_threads 4 --exon_file %s --min_samples_per_intron %d --min_samples_per_group %d %s %s --output_prefix %s"%(
            os.path.join(self.g_pipePath,'src/SplicePrep','leafcutter_ds.R'),
            strGTF2leafcutter+"_all_exons.txt.gz",
            self.g_config['leafcutter_min_samples_per_intron'],
            self.g_config['leafcutter_min_samples_per_group'],
            junc2clust_prefix+"_perind_numers.counts.gz",
            grpF,
            leafcutterDS_prefix
        ))
        cmd.append("Rscript %s --exon_file %s %s %s %s --plot_FDR %f --max_plots %d --output %s"%(
            os.path.join(self.g_pipePath,'src/SplicePrep','leafcutter_ds_plots.R'),
            strGTF2leafcutter+"_all_exons.txt.gz",
            junc2clust_prefix+"_perind_numers.counts.gz",
            grpF,
            leafcutterDS_prefix+"_cluster_significance.txt",
            self.g_config['leafcutter_plots_fdr'],
            self.g_config['leafcutter_max_plots'],
            leafcutterDS_plot
        ))
        cmd.append("Rscript %s --meta_data_file %s --output %s --FDR %f %s %s %s %s"%(
            os.path.join(self.g_pipePath,'src/SplicePrep','leafcutter_prepare_results.R'),
            grpF,
            leafcutterDS_rdata,
            self.g_config['leafviz_fdr_threshold'],
            junc2clust_prefix+"_perind_numers.counts.gz",
            leafcutterDS_prefix+"_cluster_significance.txt",
            leafcutterDS_prefix+"_effect_sizes.txt",
            strGTF2leafcutter
        ))
        return self.getOneSbatch('LC_'+strCom,strSbatch,cmd,4)

    def rmats_run(self):
        strOut,strTmp,strSbatch = self.initRun('rmats')
        if strOut is None:
            return
        for i in range(self.g_cInfo.shape[0]):
            for g in ['Alt','Ref']:
                strJob = re.sub(" ","_",self.g_cInfo.Group[i]+"_"+self.g_cInfo[g][i])
                if os.path.isfile(os.path.join(strTmp,strJob+'.txt')):
                    continue
                with open(os.path.join(strTmp,strJob+'.txt'),"w") as f:
                    f.write(','.join([os.path.join(self.g_config['bam_path'],_+"."+self.g_config['bam_suffix']) for _ in self.g_sInfo[self.g_sInfo[self.g_cInfo.Group[i]]==self.g_cInfo[g][i]]['Sample_Name']]))
        var_read_length = " --variable-read-length" if self.g_config['rmats_variable_read_length'] else ""
        paired_stats = " --paired-stats" if self.g_config['rmats_paired_stats'] else ""
        clipping = " --allow-clipping" if self.g_config['rmats_allow_clipping'] else ""
        mil = "" if self.g_config['rmats_min_intron_length_for_nss'] is None else " --mil %d"%self.g_config['rmats_min_intron_length_for_nss']
        mel = "" if self.g_config['rmats_max_exon_length_for_nss'] is None else " --mel %d"%self.g_config['rmats_max_exon_length_for_nss']
        novel_splice_sites = mil+mel if self.g_config['rmats_novel_splice_sites'] else ""
        library_type = "fr-unstranded" if self.g_config['strand']==0 else ("fr-firststrand" if self.g_config['strand']==1 else "fr-secondstrand")
        allJobs=[]
        for i in range(self.g_cInfo.shape[0]):
            oneOut = os.path.join(strTmp,self.g_cInfo.Group[i]+"_g_"+self.g_cInfo.Alt[i]+"_vs_"+self.g_cInfo.Ref[i])
            os.makedirs (oneOut,exist_ok=True)
            if len(glob.glob(os.path.join(oneOut,"*MATS.JCEC.txt")))>0:
                continue
            strRef = re.sub(" ","_",self.g_cInfo.Group[i]+"_"+self.g_cInfo.Ref[i])+".txt"
            strAlt = re.sub(" ","_",self.g_cInfo.Group[i]+"_"+self.g_cInfo.Alt[i])+".txt"
            cmd = ["module load rMATS-turbo"]
            if self.g_config['rmats_novel_splice_sites']:
                cmd.append("rmats.py --novelSS --b1 %s --b2 %s --gtf %s -t %s --libType %s --readLength %d --nthread 4 --od %s --tmp %s --cstat %f%s"%(
                    os.path.join(strTmp,strRef),
                    os.path.join(strTmp,strAlt),
                    self.g_config['genome_gtf'],
                    self.g_config['reads_type'],
                    library_type,
                    self.g_config['read_length'],
                    oneOut,
                    os.path.join(oneOut,'rmat_tmp'),
                    self.g_config['rmats_cstat'],
                    var_read_length+paired_stats+clipping+novel_splice_sites
                ))
            else:
                cmd.append("rmats.py --b1 %s --b2 %s --gtf %s -t %s --libType %s --readLength %d --nthread 4 --od %s --tmp %s --cstat %f%s"%(
                    os.path.join(strTmp,strRef),
                    os.path.join(strTmp,strAlt),
                    self.g_config['genome_gtf'],
                    self.g_config['reads_type'],
                    library_type,
                    self.g_config['read_length'],
                    oneOut,
                    os.path.join(oneOut,'rmat_tmp'),
                    self.g_config['rmats_cstat'],
                    var_read_length+paired_stats+clipping+novel_splice_sites
                ))
            cmd.append("cd %s"%oneOut)
            cmd.append("find -type f -name '*MATS.JCEC.txt' -exec sh -c 'cp $0 %s_vs_%s_$(basename $0)' {} \;"%(self.g_cInfo.Alt[i],self.g_cInfo.Ref[i]))
            oneScript = self.getOneSbatch('rM_'+os.path.basename(oneOut),strSbatch,cmd,4)
            if self.g_debug:
                print('\t%s'%os.path.basename(os.path.basename(oneScript)))
            else:
                allJobs.append(subprocess.run("sbatch %s"%oneScript,
                    shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        return allJobs

    def majiq_run(self, index_jobID):
        strOut,strTmp,strSbatch = self.initRun('majiq')
        if strOut is None:
            return
        # gtf to exon
        strGFF3 = os.path.join(strTmp,'annotation.gff3')
        jobIDs = index_jobID or []
        if not os.path.isfile(strGFF3):
            cmd = ["module load AGAT"]
            cmd.append('agat_convert_sp_gxf2gxf.pl -g %s -o %s'%(self.g_config['genome_gtf'],strGFF3))
            oneBatch = self.getOneSbatch('MJ_GFF3',strSbatch,cmd,CoreNum=2)
            if self.g_debug:
                print("\tGFF3")
            else:
                jobIDs.append(subprocess.run("sbatch %s"%oneBatch,shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        # build
        allJobs=[]
        job_dependency1 = "--dependency=afterok:%s"%','.join(jobIDs) if len(jobIDs)>0 else ""
        for row in self.g_cInfo.itertuples():
            oneBatch,oneOut = self.majiq_One(row,strGFF3,strTmp,strSbatch)
            if oneBatch is not None:
                if self.g_debug:
                    print("\t%s"%os.path.basename(oneBatch))
                else:
                    oneJobID=subprocess.run("sbatch %s %s"%(job_dependency1,oneBatch),
                        shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip()
                    job_dependency2 = "--dependency=afterok:%s"%oneJobID
            for cutoff in self.g_config['majiq_cutoffs']:
                oneScript = self.majiq_voila(oneOut,cutoff,strSbatch)
                if oneScript is None:
                    continue
                if self.g_debug:
                    print("\t%s"%os.path.basename(oneScript))
                else:
                    allJobs.append(subprocess.run("sbatch %s %s"%(job_dependency2,oneScript),
                        shell=True,check=True,capture_output=True,text=True).stdout.split(" ")[-1].strip())
        return allJobs
    def majiq_One(self, one,strGFF3,strTmp,strSbatch):
        #build
        strCom = one.Group +"_g_"+ one.Alt+"_vs_"+one.Ref
        strOut = os.path.join(strTmp,strCom)
        os.makedirs(strOut,exist_ok=True)
        if os.path.isfile(os.path.join(strOut,'reference-alternative.deltapsi.voila')):
            return None,strOut
        refBam = [_+"."+re.sub('.bam$','',self.g_config['bam_suffix']) for _ in self.g_sInfo[self.g_sInfo[one.Group]==one.Ref].Sample_Name]
        altBam = [_+"."+re.sub('.bam$','',self.g_config['bam_suffix']) for _ in self.g_sInfo[self.g_sInfo[one.Group]==one.Alt].Sample_Name]
        config = ['[info]']
        config.append('bamdirs=%s'%self.g_config['bam_path'])
        config.append('genome=%s'%self.g_config['genome_fa'])
        config.append('strandedness=forward')
        config.append('[experiments]')
        config.append('reference=%s'%','.join(refBam))
        config.append('alternative=%s'%','.join(altBam))
        strConfig = os.path.join(strOut,'majiq_v2_config.txt')
        with open(strConfig,'w') as f:
            f.write('\n'.join(config))
        cmd=["module load MAJIQ/2.2.0"]
        cmd.append('majiq build %s -c %s -o %s -j 8'%(strGFF3,strConfig,strOut))
        #deltapsi
        refMajiq = [os.path.join(strOut,_+"."+re.sub('bam$','majiq',self.g_config['bam_suffix'])) for _ in self.g_sInfo[self.g_sInfo[one.Group]==one.Ref].Sample_Name]
        altMajiq = [os.path.join(strOut,_+"."+re.sub('bam$','majiq',self.g_config['bam_suffix'])) for _ in self.g_sInfo[self.g_sInfo[one.Group]==one.Alt].Sample_Name]
        ref_str = ' '.join(refMajiq)
        alt_str = ' '.join(altMajiq)
        delta_cmd = (
            'majiq deltapsi -grp1 {ref} -grp2 {alt} -o {out} -n reference alternative -j 8 || '
            'majiq deltapsi --default-prior -grp1 {ref} -grp2 {alt} -o {out} -n reference alternative -j 8'
            ).format(ref=ref_str, alt=alt_str, out=strOut)
        cmd.append(delta_cmd)
        return self.getOneSbatch('MJ_%s'%strCom,strSbatch,cmd,CoreNum=8),strOut
    
    def majiq_voila(self, strOut,cutoff,strSbatch):
        splicegraph=os.path.join(strOut,'splicegraph.sql')
        deltapsi=os.path.join(strOut,'reference-alternative.deltapsi.voila')
        strVoila = os.path.join(strOut,os.path.basename(strOut).split("_g_")[1]+'_cutoff_%.1f_reference-alternative.deltapsi.voila.tsv'%cutoff)
        if os.path.isfile(strVoila):
            return None
        cmd=["module load MAJIQ/2.2.0",'export TMPDIR=/tmp']
        cmd.append("voila tsv %s %s -f %s --changing-between-group-dpsi %.2f --threshold %.2f --show-all -j 8"%(splicegraph,deltapsi,strVoila,cutoff,cutoff))
        return self.getOneSbatch('MJ_%s_%.2f'%(os.path.basename(strOut),cutoff),strSbatch,cmd,CoreNum=8)
    

    def wait_for_jobs_to_complete(self, timeout=14400, poll_interval=10):
        """
        Polls the log files until all submitted jobs have 'DONE' in the last line
        or until a timeout is reached.
        """
        start_time = time.time()
        while True:
            all_done = True
            for job in self.g_config['run']:
                log_files = glob.glob(os.path.join(self.g_config['output_path'], f'{job}_output', 'sbatch', '*log'))
                for log_file in log_files:
                    # Only check logs that have a corresponding sh file (job was executed)
                    sh_file = re.sub('log$', 'sh', log_file)
                    if not os.path.isfile(sh_file):
                        continue
                    with open(log_file, 'r') as f:
                        lines = f.readlines()
                        # If file is empty or doesn't end with 'DONE', mark as not done
                        if not lines or not lines[-1].startswith('DONE'):
                            all_done = False
                            break  # break inner loop to poll again
                if not all_done:
                    break  # break outer loop as well
            if all_done:
                print("All jobs are complete.")
                break
            elif time.time() - start_time > timeout:
                print("Timeout reached while waiting for jobs to complete.")
                break
            else:
                print("Waiting for jobs to complete...")
                time.sleep(poll_interval)

    
    def finalize(self, strConfig, jIDs=None):
        if jIDs is not None:
            print("Finalize job")
            if len(jIDs)>0:
                cmd=['module load Anaconda3 &> /dev/null']
                cmd.append('umask 002')
                cmd.append('python %s/src/SplicePrep.py final %s'%(self.g_pipePath,os.path.abspath(strConfig)))
                oneScript = self.getOneSbatch('SP_final',self.g_config['output_path'],cmd,1)
                if self.g_debug:
                    print("\t%s"%os.path.basename(oneScript))
                else:
                    dependency = ','.join(jIDs)
                    subprocess.run("sbatch --dependency=afterok:%s %s" % (dependency, oneScript),
                               shell=True, check=True, capture_output=True, text=True)
            # return
        
        self.getSysConfig()
        self.wait_for_jobs_to_complete()
        
        # check the submitted jobs
        print("Checking the completed jobs:")
        for job in self.g_config['run']:
            print("\t%s"%job)
            log_files = glob.glob(os.path.join(self.g_config['output_path'], f'{job}_output', 'sbatch', '*log'))
            for log_file in log_files:
                sh_file = re.sub('log$', 'sh', log_file)
                if not os.path.isfile(sh_file):
                    continue
                with open(log_file, 'r') as f:
                    lines = f.readlines()
                    if not lines[-1].startswith('DONE'):
                        print("\t\tError in", os.path.basename(log_file))
        # create projects
        print("Finalizing files ...")
        strPrj = self.g_config['wdir']
        if self.g_config['TST'] is not None:
            # timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
            # user = os.getenv('USER')
            # strPrj = os.path.join(self.g_sysConfig['localPath'], self.g_config['TST'], f"{timestamp}_{user}")
            os.makedirs(strPrj, exist_ok=True)
            shutil.copy(self.g_config['sample_info'], strPrj)
            shutil.copy(self.g_config['compare_info'], strPrj)
            self.g_config['sample_info'] = os.path.join(strPrj, os.path.basename(self.g_config['sample_info']))
            self.g_config['compare_info'] = os.path.join(strPrj, os.path.basename(self.g_config['compare_info']))
            with open(os.path.join(strPrj, 'config.yml'), "w") as f:
                yaml.dump(self.g_config, f)
                
        for job in self.g_config['run']:
            if job in self.g_sysConfig.keys():
                strLocal = os.path.join(self.g_config['output_path'],'%s_output'%job)
                onePrj = None
                if strPrj is not None:
                    onePrj=os.path.join(strPrj,job)
                    os.makedirs(onePrj,exist_ok=True)
                for pattern in self.g_sysConfig[job]:
                    pattern_path = os.path.join(strLocal, 'tmp', '**', pattern)
                    for f in glob.glob(pattern_path, recursive=True):
                        dest_local = os.path.join(strLocal, os.path.basename(f))
                        if not os.path.isfile(dest_local):
                            shutil.copy(f, strLocal)
                        if onePrj is not None:
                            dest_project = os.path.join(onePrj, os.path.basename(f))
                            if not os.path.isfile(dest_project):
                                shutil.copy(f, dest_project)
                        # print(f, dest_local, dest_project)
        if strPrj is not None:
            print("\n*** A project folder is created in: %s***\n"%strPrj)
        
        return 

    def run(self):
        # Load config from command-line argument (sys.argv[2])
        self.checkConfig()
        self.getSbatch()
        index_jobID = self.getBAMindex()
        allJobs = []
        allJobs += self.stringtie_run() or []
        allJobs += self.leafcutter_run(index_jobID) or []
        allJobs += self.rmats_run() or []
        allJobs += self.majiq_run(index_jobID) or []
        print('allJobs', allJobs)
        self.finalize(allJobs)
        return 

def main(config=None, config_file=None, strPrj=None):
    if config is not None:
        pipeline = SplicePrepPipeline(config)
        pipeline.run()
        return
    # Otherwise, use sys.argv to determine mode and config.
    mode = sys.argv[1] 
    config_file = sys.argv[2]
    config = SplicePrepPipeline.getConfig(config_file)
    pipeline = SplicePrepPipeline(config)
    if mode == 'main':
        pipeline.run()
    elif mode == 'final':
        pipeline.finalize(config_file) 
    else:
        pipeline.msgError(["Error: unknown task: %s" % mode])

if __name__ == '__main__':
    main()

