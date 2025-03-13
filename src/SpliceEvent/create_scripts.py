import argparse
import os
import shutil
import glob
import ast
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-indir', dest="indir", help ='data directory')
parser.add_argument('-bamindir', dest='bamindir', help='bam data directory', default=None)
parser.add_argument('-wdir', dest="wdir", help='working directory')
parser.add_argument('-ref', dest='ref', help="reference .gtf")
#parser.add_argument('-comparisonRef_name', dest = "comparisonRef_name", help = "reference name in comparisons")
parser.add_argument('-specific_comparison', dest= 'specific_comparison', help = "identify your specific comparison substring...", default=None)
parser.add_argument('-junction_filter', dest='junction_filter', help="open or close this filter", default=False)
parser.add_argument('-majiq_cutoff_val', dest='majiq_cutoff_val', type=float, default=0.95)
parser.add_argument('-junctionFC', dest='junctionFC', type=float, default=1.2)
parser.add_argument('-junctionMAX', dest='junctionMAX', type=float, default=25)
parser.add_argument('-samplesheet', dest='samplesheet', help= 'sample sheet file', default=None)
parser.add_argument('-addexons', dest='addexons', help= 'added exons sheet file', default='/home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/data/STMN2_add_exon.csv')

args = parser.parse_args()

args.junction_filter = args.junction_filter.lower() == 'true'

data_folder = args.indir
bam_folder =args.bamindir
wd_folder = args.wdir
reference = args.ref
specific_comparison_str = args.specific_comparison
print(specific_comparison_str)
sample_sheet_file = args.samplesheet

print(args.junction_filter)

specific_comparison = []

if specific_comparison_str:
    try:
        specific_comparison = ast.literal_eval(specific_comparison_str)
        if not isinstance(specific_comparison, list):
            raise ValueError("specific_comparison should be a list")
    except ValueError as e:
        print(f"Error in parsing specific_comparison: {e}")
        specific_comparison = []  # Reset to empty list on error

stringtie_path = os.path.join(data_folder, "stringtie")
try:
    lst_ = os.listdir(stringtie_path)
except FileNotFoundError:
    print(f"Directory not found: {stringtie_path}")
    lst_ = []

lst_NEW = [ele for ele in lst_ if any(sub_ in ele for sub_ in specific_comparison)] if specific_comparison else lst_

for i, file in enumerate(lst_NEW):
    file_prefix = file.split('/')[-1].replace('.combined.gtf', '')
    os.makedirs(wd_folder + '/' + file_prefix, exist_ok=True)
    os.makedirs(wd_folder + '/' +file_prefix + '/harm', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/harm/mxe', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/harm/mxe/out', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/harm/out', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/junction_prep', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/out', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/out/mxe', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/out/all_gene', exist_ok=True)
    os.makedirs(wd_folder + '/' + file_prefix + '/out/all_gene_junction', exist_ok=True)

    os.chdir(wd_folder + '/' + file_prefix)
    comparisonRef_name = file_prefix.split('_vs_')[-1]
    rmats_input = data_folder + '/rmats/' + file_prefix  + '_*JCEC.txt' 
    leafcutter_input= data_folder + '/leafcutter/' + file_prefix + '*_cluster_significance.txt'
    majiq_input = data_folder + '/majiq/' + file_prefix + '*voila.tsv'
    rmats_output= wd_folder  + '/' + file_prefix + '/junction_prep/rmats_junction_prep.csv' 
    leafcutter_output= wd_folder  + '/' + file_prefix + '/junction_prep/leafcutter_junction_prep.csv' 
    majiq_output= wd_folder  + '/' + file_prefix + '/junction_prep/majiq_junction_prep.csv'
    
    junction_input = wd_folder  + '/' + file_prefix + '/junction_prep/{method}_junction_prep.csv'
    sj_output = wd_folder  + '/' + file_prefix + '/junction_prep/sj.csv' 
    
    annotation= data_folder + "/stringtie/" + file
    harm_prep_dir= wd_folder  + '/' +file_prefix + "/harm"
    mxe_prep_dir = wd_folder + '/' +file_prefix + "/harm/mxe"
    harm_post_input =  harm_prep_dir + '/harm*.csv'
    mxe_post_input = mxe_prep_dir + '/harm*.csv'
    harm_post_dir= wd_folder + '/' +file_prefix + "/harm/out"
    mxe_post_dir= wd_folder + '/' +file_prefix + "/harm/mxe/out"
    relabel = harm_prep_dir + '/exon_relabel.csv'
    majiq_cutoff = wd_folder + '/' + file_prefix + '/harm/majiq_cutoff_value.csv'
    

    final_output= wd_folder  + '/' + file_prefix + "/out"
    final_mxe_output= wd_folder  + '/' + file_prefix + "/out/mxe"
    final_allgene_output=  wd_folder  + '/' + file_prefix + "/out/all_gene"
    all_gene_input =  final_allgene_output + '/*.csv'
    all_junction_output =  wd_folder  + '/' + file_prefix + "/out/all_gene_junction"
    
    bam_dir = bam_folder

    bam_pattern = "Aligned.sortedByCoord.out.bam"

    if args.junction_filter and sample_sheet_file:
        sample_sheet = pd.read_csv(sample_sheet_file   ,sep="\t")
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
        altbam_files = [bam_dir + '/' + x for x in alt_bams.to_list()]

        refbam_files = [bam_dir + '/' + x for x in ref_bams.to_list()]



    # alt_bam_str = bam_dir + '/' + file_prefix.split('_vs_')[0] + '*.Aligned.sortedByCoord.out.bam'
    # ref_bam_str = bam_dir + '/' + file_prefix.split('_vs_')[-1] + '*.Aligned.sortedByCoord.out.bam'
    # altbam_files = glob.glob(bam_dir + '/' + alt_bam_str)
    # refbam_files = glob.glob(bam_dir + '/' + ref_bam_str)
    

    inputjson=wd_folder +  '/' + file_prefix +'/out/libdepth_counts.json'

    if args.addexons:

        file_content= ['#!/bin/bash\n', 
                    '#SBATCH --job-name=harm_run_1\n',
                    '#SBATCH --nodes=1\n', 
                    '#SBATCH --cpus-per-task=12\n',
                    '#SBATCH --partition=cpu\n',
                    '#SBATCH --time=24:00:00\n',
                    '#SBATCH --mem=192G\n',
                    '#SBATCH -o harm_run_1.o.log\n',
                    '#SBATCH -e harm_run_1.e.log\n',
                    'source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh\n',
                    'conda activate splicing\n',
                    'rmats_input="{}"\n',
                    'leafcutter_input="{}"\n',
                    'majiq_input="{}"\n',
                    'rmats_output="{}"\n',
                    'leafcutter_output="{}"\n',
                    'majiq_output="{}"\n',
                    'comparisonRef_name="{}"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/rmats_junction_prep_V3.R -i "$rmats_input" -o "$rmats_output" -n "$comparisonRef_name"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/leafcutter_junction_prep_4simulation.R -i "$leafcutter_input" -o "$leafcutter_output"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/majiq_junction_prep.R -i "$majiq_input" -o "$majiq_output"\n',
                    'annotation="{}"\n',
                    'reference="{}"\n',
                    'junction_input="{}"\n',
                    'sj_output="{}"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/annotation_SJ_v3.R -a "$annotation" -r "$reference" -i "$junction_input" -o "$sj_output"\n',
                    'harm_prep_dir="{}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/harm_prep_V4_remove_mapping_V3_add_exons.py -a "$annotation" -r "$reference" -sj "$sj_output" -o "$harm_prep_dir" -addexons "{}"\n',
                    '',
                    ]
        file_content= ''.join(file_content).format(rmats_input, leafcutter_input, majiq_input, rmats_output,
                                                    leafcutter_output, majiq_output, comparisonRef_name, 
                                                    annotation, reference, junction_input, sj_output, harm_prep_dir, args.addexons)
    else:

        file_content= ['#!/bin/bash\n', 
                    '#SBATCH --job-name=harm_run_1\n',
                    '#SBATCH --nodes=1\n', 
                    '#SBATCH --cpus-per-task=12\n',
                    '#SBATCH --partition=cpu\n',
                    '#SBATCH --time=24:00:00\n',
                    '#SBATCH --mem=192G\n',
                    '#SBATCH -o harm_run_1.o.log\n',
                    '#SBATCH -e harm_run_1.e.log\n',
                    'source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh\n',
                    'conda activate splicing\n',
                    'rmats_input="{}"\n',
                    'leafcutter_input="{}"\n',
                    'majiq_input="{}"\n',
                    'rmats_output="{}"\n',
                    'leafcutter_output="{}"\n',
                    'majiq_output="{}"\n',
                    'comparisonRef_name="{}"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/rmats_junction_prep_V3.R -i "$rmats_input" -o "$rmats_output" -n "$comparisonRef_name"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/leafcutter_junction_prep_V2.R -i "$leafcutter_input" -o "$leafcutter_output"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/majiq_junction_prep.R -i "$majiq_input" -o "$majiq_output"\n',
                    'annotation="{}"\n',
                    'reference="{}"\n',
                    'junction_input="{}"\n',
                    'sj_output="{}"\n',
                    'Rscript /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/annotation_SJ_v3.R -a "$annotation" -r "$reference" -i "$junction_input" -o "$sj_output"\n',
                    'harm_prep_dir="{}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/harm_prep_V4_remove_mapping_V3_add_exons.py -a "$annotation" -r "$reference" -sj "$sj_output" -o "$harm_prep_dir"\n',
                    '',
                    ]
        file_content= ''.join(file_content).format(rmats_input, leafcutter_input, majiq_input, rmats_output,
                                                    leafcutter_output, majiq_output, comparisonRef_name, 
                                                    annotation, reference, junction_input, sj_output, harm_prep_dir)
    


    with open("run1.sh", 'w') as file:
        file.write(file_content)
        
        
    file_content = ['#!/bin/bash\n',
                    '#SBATCH --job-name=harm_run_2\n',
                    '#SBATCH --nodes=1\n',
                    '#SBATCH --partition=cpu\n',
                    '#SBATCH --time=48:00:00\n',
                    '#SBATCH --cpus-per-task=1\n',
                    '#SBATCH --mem=8GB\n',
                    '#SBATCH -o harm_run_2.o.log\n',
                    '#SBATCH -e harm_run_2.e.log\n',
                    '#SBATCH --array=1-80\n',
                    'source /home/ychen12/anaconda3/etc/profile.d/conda.sh\n',
                    'conda activate splicing\n',
                    'relabel="{}"\n',
                    'harm_post_dir="{}"\n',
                    'majiq_cutoff_file="{}"\n',
                    'input_files=($(ls {}))\n',
                    'mxe_files=($(ls {}))\n',
                    'combined_files=("${{input_files[@]}}" "${{mxe_files[@]}}")\n',
                    'input_csv="${{combined_files[$SLURM_ARRAY_TASK_ID - 1]}}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/harm_assign_job_array_V3.py -exonlabel $relabel -infile $input_csv -outdir $harm_post_dir -majiqcutofffile $majiq_cutoff_file\n',
                    ]
    file_content = ''.join(file_content).format(relabel, harm_post_dir, majiq_cutoff, harm_post_input, mxe_post_input)

    with open("run2.sh", 'w') as file:
        file.write(file_content)
    
    if args.junction_filter:
        file_content = ['#!/bin/bash\n', 
                    '#SBATCH --job-name=harm_run_3\n',
                    '#SBATCH --nodes=1\n', 
                    '#SBATCH --partition=cpu\n',
                    '#SBATCH --time=4:00:00\n',
                    '#SBATCH -o harm_run_3.o.log\n',
                    '#SBATCH -e harm_run_3.e.log\n',
                    'source /home/ychen12/anaconda3/etc/profile.d/conda.sh\n',
                    'conda activate splicing\n',
                    'file_outdir="{}"\n',
                    'final_output="{}"\n',
                    'final_allgene_output="{}"\n',
                    'mxe_outdir="{}"\n',
                    'mxe_output="{}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process.py -indir $file_outdir -outdir $final_output -allgene $final_allgene_output -majiq_score_cutoff "{}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process_4mxe.py -indir $mxe_outdir -outdir $mxe_output\n',
                    # 'altbam=($"{}"/$"{}")\n',
                    # 'refbam=($"{}"/$"{}")\n',
                    # 'altbam_str="[\"${altbam[@]}\"]"\n',
                    # 'altbam_str="${altbam_str// /\", \"}"\n',
                    # 'refbam_str="[\"${refbam[@]}\"]"\n', 
                    # 'refbam_str="${refbam_str// /\", \"}"\n'
                    # 'python /home/ychen12/splicing_test/splicing_harmonization/DIRS/libdepth_count.py -outdir $final_output -Altbam "$altbam_str" -Refbam "$refbam_str"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/libdepth_count.py -outdir $final_output -Altbam "{}" -Refbam "{}"\n',
                    '',
                    ]
        
        file_content= ''.join(file_content).format(harm_post_dir, final_output, final_allgene_output,  mxe_post_dir, final_mxe_output, args.majiq_cutoff_val, altbam_files, refbam_files)
        with open("run3.sh", 'w') as file:
            file.write(file_content)

    else:
        file_content = ['#!/bin/bash\n', 
                    '#SBATCH --job-name=harm_run_3\n',
                    '#SBATCH --nodes=1\n', 
                    '#SBATCH --partition=cpu\n',
                    '#SBATCH --time=4:00:00\n',
                    '#SBATCH -o harm_run_3.o.log\n',
                    '#SBATCH -e harm_run_3.e.log\n',
                    'source /home/ychen12/anaconda3/etc/profile.d/conda.sh\n',
                    'conda activate splicing\n',
                    'file_outdir="{}"\n',
                    'final_output="{}"\n',
                    'final_allgene_output="{}"\n',
                    'mxe_outdir="{}"\n',
                    'mxe_output="{}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process.py -indir $file_outdir -outdir $final_output -allgene $final_allgene_output -majiq_score_cutoff "{}"\n',
                    'python3 /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process_4mxe.py -indir $mxe_outdir -outdir $mxe_output\n',
                    '',
                    ]
        file_content= ''.join(file_content).format(harm_post_dir, final_output, final_allgene_output, mxe_post_dir, final_mxe_output, args.majiq_cutoff_val)
        with open("run3.sh", 'w') as file:
            file.write(file_content)

   
    # if args.junction_filter: 
    #     file_content = ['#!/bin/bash\n', 
    #                 '#SBATCH --job-name=harm_run_4\n',
    #                 '#SBATCH --nodes=1\n', 
    #                 '#SBATCH --array=1-50\n',
    #                 '#SBATCH --partition=cpu\n',
    #                 '#SBATCH --time=72:00:00\n',
    #                 '#SBATCH --mem-per-cpu=4G\n',
    #                 '#SBATCH -o harm_run_4.o.log\n',
    #                 '#SBATCH -e harm_run_4.e.log\n',
    #                 'source /home/ychen12/anaconda3/etc/profile.d/conda.sh\n',
    #                 'conda activate splicing\n',
    #                 'input_files=($(ls {}))\n',
    #                 'altbam=("${}".Aligned.sortedByCoord.out.bam)\n',
    #                 'refbam=("${}".Aligned.sortedByCoord.out.bam)\n',
    #                 'altbam_str="[${{altbam[@]}}]"\n',
    #                 'altbam_str="${{altbam_str// /, }}"\n',
    #                 'refbam_str="[${{refbam[@]}}]"\n',
    #                 'refbam_str="${{refbam_str// /, }}"\n',
    #                 'inputjson="{}"\n',
    #                 'file_output="{}"\n',
    #                 'files_per_task=$(( (${#input_files[@]} + SLURM_ARRAY_TASK_MAX - 1) / SLURM_ARRAY_TASK_MAX ))\n',
    #                 'start_idx=$(( (SLURM_ARRAY_TASK_ID - 1) * files_per_task ))\n',
    #                 'end_idx=$(( start_idx + files_per_task - 1 ))\n',
    #                 'if [ $end_idx -ge ${#input_files[@]} ]; then\n',
    #                 'end_idx=$((${#input_files[@]} - 1))\n',
    #                 'fi\n',
    #                 'for i in $(seq $start_idx $end_idx); do\n',
    #                 'input_csv="${input_files[$i]}"\n',
    #                 'python /home/ychen12/splicing_test/splicing_harmonization/DIRS/post_process_junction_filter_test.py -infile $input_csv -outdir $file_output -Altbam "$altbam_str" -Refbam "$refbam_str" -inputjson $inputjson -FCthreshold {} -maxcount {}\n',
    #                 'done\n',
    #                 '',
    #                 ]
    #     file_content= ''.join(file_content).format(all_gene_input, alt_bam_str, ref_bam_str, inputjson, all_junction_output, args.junctionFC, args.junctionMAX)
    #     with open("run4.sh", 'w') as file:
    #         file.write(file_content)

    if args.junction_filter:
        file_content = [
            '#!/bin/bash\n',
            '#SBATCH --job-name=harm_run_4\n',
            '#SBATCH --nodes=1\n',
            '#SBATCH --array=1-50\n',
            '#SBATCH --partition=cpu\n',
            '#SBATCH --time=72:00:00\n',
            '#SBATCH --mem-per-cpu=8G\n',
            '#SBATCH -o harm_run_4.o.log\n',
            '#SBATCH -e harm_run_4.e.log\n',
            'source /home/ychen12/anaconda3/etc/profile.d/conda.sh\n',
            'conda activate splicing\n',
            f'input_files=($(ls {all_gene_input}))\n',
            f'mxe_csv="{final_mxe_output}/all_cleaned.csv"\n',
            'combined_files=("${input_files[@]}" "${mxe_csv[@]}")\n',
            f'inputjson="{inputjson}"\n',
            f'file_output="{all_junction_output}"\n',
            'files_per_task=$(( (${#combined_files[@]} + SLURM_ARRAY_TASK_MAX - 1) / SLURM_ARRAY_TASK_MAX ))\n',
            'start_idx=$(( (SLURM_ARRAY_TASK_ID - 1) * files_per_task ))\n',
            'end_idx=$(( start_idx + files_per_task - 1 ))\n',
            'if [ $end_idx -ge ${#combined_files[@]} ]; then\n',
            'end_idx=$((${#combined_files[@]} - 1))\n',
            'fi\n',
            'for i in $(seq $start_idx $end_idx); do\n',
            'input_csv="${combined_files[$i]}"\n',
            f'python /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process_junction_filter_test.py -infile $input_csv -outdir $file_output -Altbam "{altbam_files}" -Refbam "{refbam_files}" -inputjson $inputjson -FCthreshold {args.junctionFC} -maxcount {args.junctionMAX}\n',
            'done\n',
            '', 
        ]
        
        file_content = ''.join(file_content)
        with open('run4.sh', 'w') as file:
            file.write(file_content)


        file_content = [
            '#!/bin/bash\n',
            '#SBATCH --job-name=harm_run_5\n',
            '#SBATCH --nodes=1\n',
            '#SBATCH --partition=cpu\n',
            '#SBATCH --time=72:00:00\n',
            '#SBATCH --mem-per-cpu=8G\n',
            '#SBATCH -o harm_run_5.o.log\n',
            '#SBATCH -e harm_run_5.e.log\n',
            'source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh\n',
            'conda activate splicing\n',
            'input_dir="{}"\n',
            'mxe_dir="{}"\n'
            'python /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process_MXE_junction.py -indir $input_dir -mxeindir $mxe_dir -outdir $mxe_dir\n',
            '', 
        ]
        
        file_content = ''.join(file_content).format(all_junction_output, final_mxe_output)
        with open('run5.sh', 'w') as file:
            file.write(file_content)

    else:
        file_content = [
            '#!/bin/bash\n',
            '#SBATCH --job-name=harm_run_4\n',
            '#SBATCH --nodes=1\n',
            '#SBATCH --partition=cpu\n',
            '#SBATCH --time=72:00:00\n',
            '#SBATCH --mem-per-cpu=8G\n',
            '#SBATCH -o harm_run_4.o.log\n',
            '#SBATCH -e harm_run_4.e.log\n',
            'source /home/ychen12/tools/miniconda3/etc/profile.d/conda.sh\n',
            'conda activate splicing\n',
            'input_dir="{}"\n'
            'mxe_dir="{}"\n'
            'python /home/ychen12/splicing_test/edge_splicing/splicing_harmonization/DIRS/post_process_MXE.py -indir $input_dir -mxeindir $mxe_dir -outdir $mxe_dir\n',
            '', 
            ]
            
        file_content = ''.join(file_content).format(final_allgene_output, final_mxe_output)
        with open('run4.sh', 'w') as file:
            file.write(file_content)


            
    os.chdir(wd_folder)
