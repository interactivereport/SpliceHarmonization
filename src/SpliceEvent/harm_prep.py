import numpy as np
import itertools
import pandas as pd
import warnings
import re
import genomicranges as grange
import pyranges as pr
import pandas as pd
from collections import defaultdict, Counter
warnings.filterwarnings('ignore')
# from itertools import combinations
import ast
from utils import *
import argparse
from multiprocessing import Pool 
import itertools
import os
from copy import deepcopy

parser = argparse.ArgumentParser()
parser.add_argument('-r', dest="reference", help ='reference gtf')
parser.add_argument('-a', dest='annotation', help="stringtie annotation gtf")
parser.add_argument('-sj', dest="sj",  help="csv path")
parser.add_argument('-o', dest="outdir", help="file path")
parser.add_argument('-addexons', dest="addexons", help="file path",  default=None)
args = parser.parse_args()
##############
# junction library generation
ref = pr.read_gtf(args.reference)
gtr = pr.read_gtf(args.annotation)
annotation_sj = pd.read_csv(args.sj)
if args.addexons:
    add_exons = pd.read_csv(args.addexons)
    add_exon_new = pr.PyRanges(add_exons)


ref_new = ref[['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
       'Frame', 'transcript_id', 'gene_id', 'exon_number',  'gene_name', 'cmp_ref',
       'cmp_ref_gene']]
gtr_new = gtr[['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
       'Frame', 'transcript_id', 'gene_id', 'exon_number',  'gene_name', 'cmp_ref',
       'cmp_ref_gene']]


if args.addexons:
    gr_all = pr.concat([ref_new, gtr_new, add_exon_new])
else:
    gr_all = pr.concat([ref_new, gtr_new])

MaskTranscriptIDToGeneName = dict(zip(gr_all[gr_all.Feature == "transcript"].transcript_id, gr_all[gr_all.Feature == "transcript"].gene_name))
MaskTranscriptIDToGeneNameNA = { k:v for k, v in MaskTranscriptIDToGeneName.items() if str(v) == 'nan'}
MaskTranscriptIDToGeneNameAll = { k:v for k, v in MaskTranscriptIDToGeneName.items() if str(v) != 'nan'}
del gr_all

ref_new = ref_new[ref_new.transcript_id.isin(list(MaskTranscriptIDToGeneNameAll.keys()))]
gtr_new = gtr_new[gtr_new.transcript_id.isin(list(MaskTranscriptIDToGeneNameAll.keys()))]
# add_exon_new = add_exon_new[add_exon_new.transcript_id.isin(list(MaskTranscriptIDToGeneNameAll.keys()))]


gene_name_assign_df = pd.DataFrame({'gene_id': gtr_new[(gtr_new.Feature == 'transcript')][['gene_id', 'gene_name']].gene_id,
            'gene_name': gtr_new[(gtr_new.Feature == 'transcript')][['gene_id', 'gene_name']].gene_name})
gene_name_assign = gene_name_assign_df.groupby('gene_id')['gene_name'].unique().apply(list).to_dict()

ref_df = pd.DataFrame({'Chromosome':ref_new.Chromosome, 'Start':ref_new.Start, 'Feature':ref_new.Feature,
                        'End':ref_new.End, 'Strand':ref_new.Strand,
                       'transcript_id': ref_new.transcript_id, 'gene_id': ref_new.gene_id,
                       'gene_name':ref_new.gene_name, 
                      'exon_number':ref_new.exon_number})

gr_df = pd.DataFrame({'Chromosome':gtr_new.Chromosome, 'Start':gtr_new.Start, 'Feature':gtr_new.Feature,
                        'End':gtr_new.End, 'Strand':gtr_new.Strand,
                       'transcript_id': gtr_new.transcript_id, 'gene_id': gtr_new.gene_id,
                       'gene_name':gtr_new.gene_name, 
                      'exon_number':gtr_new.exon_number})
if args.addexons:
    add_exon_new_df = pd.DataFrame({'Chromosome':add_exon_new.Chromosome, 'Start':add_exon_new.Start, 'Feature':add_exon_new.Feature,
                            'End':add_exon_new.End, 'Strand':add_exon_new.Strand,
                        'transcript_id': add_exon_new.transcript_id, 'gene_id': add_exon_new.gene_id,
                        'gene_name':add_exon_new.gene_name, 
                        'exon_number':add_exon_new.exon_number})
del ref_new
del gtr_new
gr_df = gr_df[(gr_df['Feature'] == 'exon')]
ref_df = ref_df[(ref_df['Feature'] == 'exon')]

gr_df['gene_name'] = gr_df['gene_id'].map(gene_name_assign)
gr_df = gr_df.explode('gene_name')

if args.addexons:
    df_exon = pd.concat([ref_df, gr_df, add_exon_new_df])
else:
    df_exon = pd.concat([ref_df, gr_df])
    
del ref_df
del gr_df

grouped = df_exon.groupby("gene_name")
data = []
for gene, gene_ranges in grouped:
    gene_ranges = gene_ranges.sort_values(by=['Start', 'End'])
    gene_ranges["unified_exon"] = gene_ranges.groupby(['Start', 'End']).ngroup() + 1
    data.append(gene_ranges)
del df_exon

df_exon_unified = pd.concat(data)
df_exon_unified.to_csv(args.outdir + '/exon_relabel.csv', index=False)

junction_locations = {}
for gene_name, exons in  df_exon_unified.groupby(["gene_name"]):
    exons = exons.drop_duplicates('unified_exon').reset_index(drop=True)
    start_positions = exons.Start
    end_positions = exons.End
    exon_number = exons.unified_exon
    junctions = [
        tuple(z) for n in range(1, len(end_positions))
        for z in zip(end_positions[:-n] + 1, start_positions[n:],
                     [ str(x) + '-' + str(y) for x, y in zip(exon_number[n:], exon_number[:-n])])]
    # print(unified_number)
    junction_locations[gene_name] = junctions

flattened_junctions =  [(key, *value) for key, values in junction_locations.items() for value in values]
df_junction = pd.DataFrame(flattened_junctions, columns=['gene_name', 'junction_start', 'junction_end', 'unified_junction_label'])
df_junction.to_csv(args.outdir + '/permutative_junctions.csv', index=False)


del junction_locations
del flattened_junctions


########
print('junction library generation completed ...')

########
# mapping junctions from three algorithm to reference junciton library
sj = pd.read_csv(args.sj)

unique_gene_junction = {(row.gene_name, row.junction_start, row.junction_end) for row in df_junction.itertuples(index=False)}
sj['junction_filter'] = list(zip(sj['geneSymbol'], sj['start'], sj['end']))
filtered_sj = sj[sj['junction_filter'].apply(lambda x: x in unique_gene_junction)]
residual_sj = sj[~sj['junction_filter'].apply(lambda x: x in unique_gene_junction)]


residual_sj.to_csv(args.outdir + '/residual_sj.csv', index=False)
annotated_sj = residual_sj.merge(df_junction, left_on=['start', 'end', 'geneSymbol'], 
                               right_on=['junction_start', 'junction_end', 'gene_name'], 
                               how='left')

missing_gene = set(sj['geneSymbol']) - set(annotated_sj['geneSymbol'])
print('the number of missing genes by briefly mapping to junction library is {}'.format(len(missing_gene)))
print(missing_gene)

annotated_sj['exon_list'] = annotated_sj['unified_junction_label'].apply(lambda x: np.nan if pd.isnull(x) else x.split("-")[::-1])
annotated_sj['exon_list'] = annotated_sj['exon_list'].apply(lambda x: [int(item) for item in x])

rmats = annotated_sj[annotated_sj['method'] == 'rmats'].reset_index(drop=True)
leafcutter = annotated_sj[annotated_sj['method'] == 'leafcutter'].reset_index(drop=True)
majiq = annotated_sj[annotated_sj['method'] == 'majiq'].reset_index(drop=True)

# rmats.to_csv(args.outdir + '/rmats_harm_prep.csv', index=False)
# leafcutter.to_csv(args.outdir + '/leafcutter_harm_prep.csv', index=False)
# majiq.to_csv(args.outdir + '/majiq_harm_prep.csv', index=False)

# del rmats
# del leafcutter
# del majiq
del annotated_sj
del sj
del unique_gene_junction
del filtered_sj

##### rmats mapping
#rmats = pd.read_csv(args.outdir + '/rmats_harm_prep.csv')
if len(rmats)> 0:
    rmats['ID2'] = rmats.apply(lambda row: re.sub(r'^ES_SJ', 'ES', row['ID']) if row['type'] == 'ES_SJ'
                                    else re.sub(r'^EIL_SJ', 'ES', row['ID']) if row['type'] == 'EIL_SJ'
                                    else re.sub(r'^EIR_SJ', 'ES', row['ID']) if row['type'] == 'EIR_SJ'
                                    else re.sub(r'^MXE1L_SJ', 'MXE', row['ID']) if row['type'] == 'MXE1L_SJ'
                                    else re.sub(r'^MXE1R_SJ', 'MXE', row['ID']) if row['type'] == 'MXE1R_SJ'
                                    else re.sub(r'^MXE2L_SJ', 'MXE', row['ID']) if row['type'] == 'MXE2L_SJ'
                                    else re.sub(r'^MXE2R_SJ', 'MXE', row['ID']) if row['type'] == 'MXE2R_SJ'
                                    else re.sub(r'^IRL', 'IR', row['ID']) if row['type'] == 'IRL'
                                    else re.sub(r'^IRR', 'IR', row['ID']) if row['type'] == 'IRR'
                                    else re.sub(r'^A5SS_SSJ', 'A5SS', row['ID']) if row['type'] == 'A5SS_SSJ'
                                    else re.sub(r'^A5SS_ISJ', 'A5SS', row['ID']) if row['type'] == 'A5SS_ISJ'
                                    else re.sub(r'^A3SS_SSJ', 'A3SS', row['ID']) if row['type'] == 'A3SS_SSJ'
                                    else re.sub(r'^A3SS_ISJ', 'A3SS', row['ID']) if row['type'] == 'A3SS_ISJ'
                                    else row['ID'],
                            axis=1)

    rmats_ES = rmats[rmats['ID2'].str.contains('ES')]
    es_groups = rmats_ES.groupby('ID2')
    casstte_list ={}
    es_groups
    for es, sub_df in es_groups:
        #print(es)
        coverage_list = set(tuple(sorted([int(x_) for x_ in x])) for x in sub_df[sub_df['ID'].str.contains('ES_SJ')]['exon_list'].to_numpy())
        left_list = sub_df[sub_df['ID'].str.contains('EIL_SJ')]['exon_list'].to_numpy()
        right_list = sub_df[sub_df['ID'].str.contains('EIR_SJ')]['exon_list'].to_numpy()
        # es_list[es] = [coverage_list, left_list, right_list]
        listR_dict = {key: [value for _, value in group] for key, group in itertools.groupby(right_list, key=lambda x: x[0])}
        
        output = [subleft_list + [value] for subleft_list in left_list for value in listR_dict.get(subleft_list[-1], [])]
        #print(output)
        output = list(set(tuple(sorted([int(x_) for x_  in x])) for x in output))
    # print(output)
        if coverage_list:
            cleaned_output = [o_ for o_ in output if any(set(template).issubset(o_) for template in coverage_list)]
        else:
            cleaned_output = output
        
        casstte_list[es] =  [list(sublist) for sublist in cleaned_output]
    rmats_ES['graph'] = rmats_ES['ID2'].map(casstte_list)
    #rmats_ES.dropna(subset=['graph'])
    rmats_ES = rmats_ES[rmats_ES['graph'].apply(lambda x: len(x) > 0)]


    rmats_MXE = rmats[rmats['ID'].str.contains('MXE')]
    #rmats_MXE = rmats_MXE[rmats_MXE['ID2'] == 'MXE9589']
    mxe_groups = rmats_MXE.groupby('ID2')
    casstte_list ={}
    for es, sub_df in mxe_groups:
        #print(es)
        # coverage_list = set(tuple(x) for x in sub_df[sub_df['ID'].str.contains('ES_SJ')]['exon_list'].to_numpy())
        left_list_1 =  sub_df[sub_df['ID'].str.contains('MXE1L')]['exon_list'].to_numpy()
        right_list_1 = sub_df[sub_df['ID'].str.contains('MXE1R')]['exon_list'].to_numpy()
        # print(left_list_1)
        # print(right_list_1)
        #print(left_list_1.size)
        if left_list_1.size > 0 and right_list_1.size > 0:  # Check if arrays have non-zero size
            comb1 = left_right_combine(left_list_1, right_list_1)
        elif left_list_1.size > 0:
            comb1 = left_list_1
        elif right_list_1.size > 0:
            comb1 = right_list_1
        else:
            comb1 = np.array([])
        # print(comb1)    
        left_list_2 =  sub_df[sub_df['ID'].str.contains('MXE2L')]['exon_list'].to_numpy()
        right_list_2 = sub_df[sub_df['ID'].str.contains('MXE2R')]['exon_list'].to_numpy()
        # print(right_list_2)
        
        if left_list_2.size > 0 and right_list_2.size > 0:  # Check if arrays have non-zero size
            comb2 = left_right_combine(left_list_2, right_list_2)
        elif left_list_2.size > 0:
            comb2 = left_list_2
        elif right_list_2.size > 0:
            comb2 = right_list_2
        else:
            comb2 = np.array([])
        # es_list[es] = [coverage_list, left_list, right_list]
        
        # print(comb1)
        # print(comb2)
        output = set()
        if comb1 is not None and comb2 is not None:
            for subList1 in comb1:
                for subList2 in comb2:
                    # print(subList1)
                    # print(subList2)
                    merged_sublist = np.concatenate((subList1, subList2))
                    # print(merged_sublist)
                    unique_sublist = sorted([int(x) for x in np.unique(merged_sublist)])
                    #print(unique_sublist)
                    if len(unique_sublist) == 4:
                        output.add(tuple(unique_sublist))
        else:
            output = []
        # print(merged_list)
        # merged_list = list(set(merged_list))
        casstte_list[es] = [list(sublist) for sublist in output]

    rmats_MXE['graph'] = rmats_MXE['ID2'].map(casstte_list)
    #print(len(rmats_MXE))
    rmats_MXE = rmats_MXE[rmats_MXE['graph'].apply(lambda x: len(x) > 0)]


    rmats_ASS = rmats[rmats['ID'].str.contains(r'A5SS|A3SS')]
    #rmats_ASS = rmats_ASS[rmats_ASS['ID2'] == 'A3SS0']
    ASS_groups = rmats_ASS.groupby('ID2')
    casstte_list ={}
    for es, sub_df in ASS_groups:
        left_list = sub_df[sub_df['ID'].str.contains('SSJ')]['exon_list'].to_numpy()
        right_list = sub_df[sub_df['ID'].str.contains('ISJ')]['exon_list'].to_numpy()
        # es_list[es] = [coverage_list, left_list, right_list]
        output = []
        for left_ in left_list:
            for right_ in right_list:
                if left_[0] == right_[0] or left_[-1] == right_[-1]:
                    temp_ = left_ + right_
                    output.append(list(set(temp_))) 
        output = list(set(tuple(sorted([int(x_) for x_  in x])) for x in output))
        casstte_list[es] =  [list(sublist) for sublist in output]
    rmats_ASS['graph'] = rmats_ASS['ID2'].map(casstte_list)
    #print(len(rmats_ASS))
    rmats_ASS = rmats_ASS[rmats_ASS['graph'].apply(lambda x: len(x) > 0)]

    rmats_IR = rmats[rmats['ID2'].str.contains('IR')]
    rmats_IR['graph'] = rmats_IR['exon_list'].apply(lambda x: [[int(x_) for x_ in x]])
    redefine_rmats_type = pd.concat([rmats_ES, rmats_MXE, rmats_ASS, rmats_IR])

    psi1 = redefine_rmats_type['PSI_1'].apply(lambda x: [float(val) for val in x.split(',') if val != 'NA'])
    redefine_rmats_type['PSI_1'] = [np.mean(val) for val in psi1]
    psi2 = redefine_rmats_type['PSI_2'].apply(lambda x: [float(val) for val in x.split(',') if val != 'NA'])
    redefine_rmats_type['PSI_2'] = [np.mean(val) for val in psi2]

    df_repeated = pd.DataFrame(np.repeat(redefine_rmats_type.values, [len(x) for x in redefine_rmats_type['graph']], axis=0), columns=redefine_rmats_type.columns)
    nested_cols = pd.DataFrame(redefine_rmats_type['graph'].explode().reset_index(drop=True))
    df_repeated['graph'] = nested_cols
    df_repeated['casstte_list'] = df_repeated.apply(lambda row: find_triplets_graph(row['graph']), axis=1)

    df_repeated_casstte = pd.DataFrame(np.repeat(df_repeated.values, [len(x) for x in df_repeated['casstte_list']], axis=0), columns=df_repeated.columns)
    nested_cols_casstte = pd.DataFrame(df_repeated['casstte_list'].explode().reset_index(drop=True))
    df_repeated_casstte['casstte_list'] = nested_cols_casstte

    # df_repeated = pd.DataFrame(np.repeat(redefine_rmats_type.values, [len(x) for x in redefine_rmats_type['graph']], axis=0), columns=redefine_rmats_type.columns)
    # nested_cols = pd.DataFrame(redefine_rmats_type['graph'].explode().reset_index(drop=True))
    # df_repeated['graph'] = nested_cols
    # df_repeated['casstte_list'] = df_repeated.apply(lambda row: find_triplets_graph(row['graph']), axis=1)

    # df_repeated_casstte = pd.DataFrame(np.repeat(df_repeated.values, [len(x) for x in df_repeated['casstte_list']], axis=0), columns=df_repeated.columns)
    # nested_cols_casstte = pd.DataFrame(df_repeated['casstte_list'].explode().reset_index(drop=True))
    # df_repeated_casstte['casstte_list'] = nested_cols_casstte

    casstte_array = df_repeated_casstte['casstte_list'].tolist()
    exon_array = df_repeated_casstte['exon_list'].tolist()
    is_sublist = [np.all(np.isin(sublist, mainlist)) for sublist, mainlist in zip(exon_array, casstte_array)]

    df_repeated_casstte['is_sublist'] = is_sublist
    df_repeated_casstte_cleaned = df_repeated_casstte[df_repeated_casstte['is_sublist'] == True ].reset_index(drop=True)

    df_repeated_casstte_cleaned = df_repeated_casstte_cleaned.reset_index(drop=True)
    df_repeated_casstte_cleaned['ID'] = df_repeated_casstte_cleaned['ID2'] 
    df_repeated_casstte_cleaned = df_repeated_casstte_cleaned.drop(columns=['ID2'])
    df_repeated_casstte_cleaned['exon_list'] = df_repeated_casstte_cleaned['exon_list'].apply(str)
    df_repeated_casstte_cleaned['graph'] = df_repeated_casstte_cleaned['graph'].apply(str)
    # df_repeated_casstte_cleaned = df_repeated_casstte_cleaned.drop_duplicates(subset=['ID', 'coordinates', 'score', 'comparison', 'method', 'PSI_1', 'PSI_2',
    #        'gene_id',  'gene_name2', 'exon_list'])
    df_repeated_casstte_cleaned['exon_list'] = df_repeated_casstte_cleaned['exon_list'].apply(ast.literal_eval)
    df_repeated_casstte_cleaned['graph'] = df_repeated_casstte_cleaned['graph'].apply(ast.literal_eval)

    df_repeated_casstte_cleaned['split_id'] = df_repeated_casstte_cleaned['ID'] 
    non_mxe = df_repeated_casstte_cleaned[~ df_repeated_casstte_cleaned['ID'].str.contains('MXE')]
    mxe = df_repeated_casstte_cleaned[df_repeated_casstte_cleaned['ID'].str.contains('MXE')]

    split_id = non_mxe['split_id'].unique()
    split_parts = np.array_split(split_id, 10)

    rmats_part_index = 1
    for sp in split_parts:
        sub_df = non_mxe[non_mxe['split_id'].isin(sp)]
        sub_df.to_csv(args.outdir + '/harm_rmats_p' + str(rmats_part_index) + '.csv')
        rmats_part_index += 1

    split_id = mxe['split_id'].unique()
    split_parts = np.array_split(split_id, 10)
    rmats_part_index = 1
    for sp in split_parts:
        sub_df = mxe[mxe['split_id'].isin(sp)]
        mxe_dir = args.outdir + '/mxe'
        sub_df.to_csv(mxe_dir  + '/harm_rmats_p' + str(rmats_part_index) + '.csv')
        rmats_part_index += 1

    print('rmats junction mapping done ....')
    missing_gene = set(rmats['gene_name']) - set(df_repeated_casstte_cleaned['gene_name'])
    print('the number of rmats total conatins {}'.format(len(set(rmats['gene_name']))))
    print('the number of missing genes by rmats mapping junction library is {}'.format(len(missing_gene)))
    print(missing_gene)

    del rmats
    del df_repeated_casstte
    del df_repeated_casstte_cleaned
    del non_mxe
    del mxe
    del df_repeated

##### leafcutter mapping
# leafcutter = pd.read_csv(args.outdir + '/leafcutter_harm_prep.csv')
if len(leafcutter)> 0:
    redefine_leafcutter_type = leafcutter.groupby('ID').agg(list)
    redefine_leafcutter_type['graph'] =  redefine_leafcutter_type['exon_list'].apply(lambda x: sorted(list({item for sublist in x for item in sublist})))
    redefine_leafcutter_type['casstte_list'] = redefine_leafcutter_type.apply(lambda row: exonlist_to_threenode_graph(row['exon_list'], row['graph']), axis=1)
    redefine_leafcutter_type = redefine_leafcutter_type.reset_index()

    leafcutter_after_cluster = sorted(redefine_leafcutter_type['ID'].apply(lambda x: int(x.split('_')[-2])).unique())
    leafcutter_before_cluster = sorted(leafcutter['ID'].apply(lambda x: int(x.split('_')[-2])).unique())

    leafcutter_repeat = redefine_leafcutter_type.explode(['score', 'comparison', 'method', 'PSI_1', 'PSI_2',
        'type', 'geneSymbol', 'start', 'end', 'chr', 'strand',
        'coordinates', 'junction_filter', 'gene_name', 'junction_start',
        'junction_end', 'unified_junction_label', 'exon_list'])
    leafcutter_redefine = leafcutter_repeat.reset_index(drop=True)

    # values = leafcutter_redefine['casstte_list'].to_numpy()
    # exon_lists = leafcutter_redefine['exon_list'].to_numpy()
    # # keep_or_remove = [0]* leafcutter_redefine.shape[0]
    # # casstte_list_array = [None] * leafcutter_redefine.shape[0]
    # keep_or_remove = {}
    # casstte_list_array = {}
    # for i in range(leafcutter_redefine.shape[0]):
    #     current_graph = values[i]
    #     current_exon_list = [int(element) for element in exon_lists[i]]
    #     for edge_ in current_graph:
    #         if set(current_exon_list).issubset(edge_):
    #             current_exon_tuple = tuple(current_exon_list) 
    #             if current_exon_tuple not in keep_or_remove:
    #                 keep_or_remove[current_exon_tuple] = 1
    #             else:
    #                 keep_or_remove[current_exon_tuple] += 1
    #             if current_exon_tuple not in casstte_list_array:
    #                 casstte_list_array[current_exon_tuple] = []
    #             casstte_list_array[current_exon_tuple].append(edge_)
    # keep_idx = [i for i in range(len(keep_or_remove)) if keep_or_remove[i] == 1 ]
    # leafcutter_redefine_cleaned = leafcutter_redefine[leafcutter_redefine.index.isin(keep_idx)]
    # casstte_list_array_cleaned = [casstte_list_array[i] for i in keep_idx]
    # leafcutter_redefine_cleaned['casstte_list'] = casstte_list_array_cleaned
    def repeat_exon_list(keys):
        return [keys] * keep_or_remove.get(tuple(keys), 0) 
    def map_casstte_list(keys):
        return casstte_list_array.get(tuple(keys), 0) 


    leafcutter_all_cleaned = []
    for ID, sub_df in leafcutter_redefine.groupby('ID'):
        values = sub_df['casstte_list'].to_numpy()
        exon_lists = sub_df['exon_list'].to_numpy()
        keep_or_remove = {}
        casstte_list_array = {}
        for i in range(sub_df.shape[0]):
            current_graph = values[i]
            current_exon_list = [int(element) for element in exon_lists[i]]
            for edge_ in current_graph:
                if set(current_exon_list).issubset(edge_):
                    current_exon_tuple = tuple(current_exon_list) 
                    if current_exon_tuple not in keep_or_remove:
                        keep_or_remove[current_exon_tuple] = 1
                    else:
                        keep_or_remove[current_exon_tuple] += 1
                    if current_exon_tuple not in casstte_list_array:
                        casstte_list_array[current_exon_tuple] = []
                    casstte_list_array[current_exon_tuple].append(edge_)
        
        sub_df['repeated_exon'] = sub_df['exon_list'].apply(lambda x: repeat_exon_list(x))
        sub_df['repeated_casstte'] = sub_df['exon_list'].apply(lambda x: map_casstte_list(x))
        sub_df['repeat_exon_len'] = sub_df['repeated_exon'].apply(lambda x: len(x))
        sub_df['repeat_casstee_len'] = sub_df['repeated_casstte'].apply(lambda x: len(x) if isinstance(x, list) else 0)
        sub_df_cleaned = sub_df[sub_df['repeat_exon_len'] > 0]
        try:
            sub_df_cleaned = sub_df_cleaned.explode(['repeated_exon', 'repeated_casstte'])
        except:
            print(ID)
        leafcutter_all_cleaned.append(sub_df_cleaned)


    leafcutter_redefine_cleaned = pd.concat(leafcutter_all_cleaned)


    leafcutter_redefine_cleaned['exon_list'] = leafcutter_redefine_cleaned['repeated_exon']
    leafcutter_redefine_cleaned['casstte_list'] = leafcutter_redefine_cleaned['repeated_casstte']

    leafcutter_redefine_cleaned = leafcutter_redefine_cleaned.reset_index(drop=True)
    leafcutter_redefine_cleaned['exon_list'] = leafcutter_redefine_cleaned['exon_list'].apply(str)
    leafcutter_redefine_cleaned['graph'] = leafcutter_redefine_cleaned['graph'].apply(str)

    leafcutter_redefine_cleaned['exon_list'] = leafcutter_redefine_cleaned['exon_list'].apply(ast.literal_eval)
    leafcutter_redefine_cleaned['graph'] = leafcutter_redefine_cleaned['graph'].apply(ast.literal_eval)

    leafcutter_redefine_cleaned['split_id'] = leafcutter_redefine_cleaned['ID'] 
    split_id = leafcutter_redefine_cleaned['split_id'].unique()
    split_parts = np.array_split(split_id, 10)

    leafcutter_part_index = 1
    for sp in split_parts:
        sub_df = leafcutter_redefine_cleaned[leafcutter_redefine_cleaned['split_id'].isin(sp)]
        sub_df.to_csv( args.outdir + '/harm_leafcutter_p' + str(leafcutter_part_index) + '.csv')
        leafcutter_part_index += 1

    print('leafcutter junction mapping done ....')
    missing_gene = set(leafcutter['gene_name']) - set(leafcutter_redefine_cleaned['gene_name'])
    print('the number of leafcutter total conatins {}'.format(len(set(leafcutter['gene_name']))))
    print('the number of missing genes by leafcutter mapping junction library is {}'.format(len(missing_gene)))
    print(missing_gene)

    del leafcutter 
    del redefine_leafcutter_type
    del leafcutter_redefine_cleaned
    del leafcutter_redefine
    del leafcutter_all_cleaned


# majiq
# majiq = pd.read_csv(args.outdir + '/majiq_harm_prep.csv', index=False)
majiq['cutoff'] = majiq['comparison'].apply(lambda x: x.split('_')[-1])
majiq['cutoff'] = majiq['cutoff'].astype(float)
majiq['score'] = majiq['score'].astype(float)

majiq.to_csv(args.outdir + '/majiq_cutoff_value.csv', index=False)

cutoff_counts = majiq['cutoff'].value_counts()
if cutoff_counts.nunique() == 1:
    print("pass")
else:
    print("fail")
    raise ValueError("Error: Not all cutoff counts are equal.")

print(cutoff_counts)


# majiq_cutoff_select =  majiq.groupby('ID').filter(lambda x: (x['score'] >= 0.9).any())
# majiq = pd.concat([majiq[majiq['cutoff'] == 0.1],  majiq_cutoff_select[majiq_cutoff_select['cutoff'] > 0.1]])
# # majiq = pd.concat([majiq[majiq['cutoff'] == 0.1],  majiq.groupby('ID').filter(lambda x: (x['score'] >= 0.9).any())])
# majiq.sort_values(['cutoff', 'geneSymbol'], inplace=True, axis=0, ascending=[True, True])
majiq_default = majiq[majiq['cutoff'] == 0.1].reset_index(drop=True)
# redefine_majiq_type = majiq.groupby('ID').agg(list)
redefine_majiq_type = majiq_default.groupby('ID').agg(list)
redefine_majiq_type['graph'] =  redefine_majiq_type['exon_list'].apply(lambda x: sorted(list({item for sublist in x for item in sublist})))
redefine_majiq_type['casstte_list'] = redefine_majiq_type.apply(lambda row: exonlist_to_threenode_graph(row['exon_list'], row['graph']), axis=1)
redefine_majiq_type = redefine_majiq_type.reset_index()
majiq_repeat = redefine_majiq_type.explode(['score', 'comparison', 'method', 'PSI_1', 'PSI_2',
       'type', 'geneSymbol', 'start', 'end', 'chr', 'strand', 'cutoff',
       'coordinates', 'junction_filter', 'gene_name', 'junction_start',
       'junction_end', 'unified_junction_label', 'exon_list'])
majiq_redefine = majiq_repeat.reset_index(drop=True)

# values = majiq_redefine['casstte_list'].to_numpy()
# exon_lists = majiq_redefine['exon_list'].to_numpy()
# keep_or_remove = {}
# casstte_list_array = {}
# for i in range(majiq_redefine.shape[0]):
#     current_graph = values[i]
#     current_exon_list = [int(element) for element in exon_lists[i]]
#     for edge_ in current_graph:
#         if set(current_exon_list).issubset(edge_):
#             current_exon_tuple = tuple(current_exon_list) 
#             if current_exon_tuple not in keep_or_remove:
#                 keep_or_remove[current_exon_tuple] = 1
#             else:
#                 keep_or_remove[current_exon_tuple] += 1
#             if current_exon_tuple not in casstte_list_array:
#                 casstte_list_array[current_exon_tuple] = []
#             casstte_list_array[current_exon_tuple].append(edge_)

def repeat_exon_list(keys):
    return [keys] * keep_or_remove.get(tuple(keys), 0) 
def map_casstte_list(keys):
    return casstte_list_array.get(tuple(keys), 0) 

majiq_all_cleaned = []
for ID, sub_df in majiq_redefine.groupby(['ID', 'cutoff']):
    values = sub_df['casstte_list'].to_numpy()
    exon_lists = sub_df['exon_list'].to_numpy()
    keep_or_remove = {}
    casstte_list_array = {}
    for i in range(sub_df.shape[0]):
        current_graph = values[i]
        current_exon_list = [int(element) for element in exon_lists[i]]
        for edge_ in current_graph:
            if set(current_exon_list).issubset(edge_):
                current_exon_tuple = tuple(current_exon_list) 
                if current_exon_tuple not in keep_or_remove:
                    keep_or_remove[current_exon_tuple] = 1
                else:
                    keep_or_remove[current_exon_tuple] += 1
                if current_exon_tuple not in casstte_list_array:
                    casstte_list_array[current_exon_tuple] = []
                casstte_list_array[current_exon_tuple].append(edge_)
    
    sub_df['repeated_exon'] = sub_df['exon_list'].apply(lambda x: repeat_exon_list(x))
    sub_df['repeated_casstte'] = sub_df['exon_list'].apply(lambda x: map_casstte_list(x))
    sub_df['repeat_exon_len'] = sub_df['repeated_exon'].apply(lambda x: len(x))
    sub_df['repeat_casstee_len'] = sub_df['repeated_casstte'].apply(lambda x: len(x) if isinstance(x, list) else 0)
    sub_df_cleaned = sub_df[sub_df['repeat_exon_len'] > 0]
    try:
        sub_df_cleaned = sub_df_cleaned.explode(['repeated_exon', 'repeated_casstte'])
    except:
        print(ID)
    majiq_all_cleaned.append(sub_df_cleaned)


majiq_redefine_cleaned = pd.concat(majiq_all_cleaned)

# majiq_redefine['repeated_exon'] = majiq_redefine['exon_list'].apply(lambda x: repeat_exon_list(x))
# majiq_redefine['repeated_casstte'] = majiq_redefine['exon_list'].apply(lambda x: map_casstte_list(x))
# majiq_redefine['repeat_exon_len'] = majiq_redefine['repeated_exon'].apply(lambda x: len(x))
# majiq_redefine['repeat_casstee_len'] = majiq_redefine['repeated_casstte'].apply(lambda x: len(x) if isinstance(x, list) else 0)
# majiq_redefine_cleaned = majiq_redefine[majiq_redefine['repeat_exon_len'] > 0]
# majiq_redefine_cleaned = majiq_redefine_cleaned.explode(['repeated_exon', 'repeated_casstte'])
majiq_redefine_cleaned['exon_list'] = majiq_redefine_cleaned['repeated_exon']
majiq_redefine_cleaned['casstte_list'] = majiq_redefine_cleaned['repeated_casstte']
# keep_idx = [i for i in range(len(keep_or_remove)) if keep_or_remove[i] == 1 ]
# majiq_redefine_cleaned = majiq_redefine[majiq_redefine.index.isin(keep_idx)]
# casstte_list_array_cleaned = [casstte_list_array[i] for i in keep_idx]
# majiq_redefine_cleaned['casstte_list'] = casstte_list_array_cleaned

majiq_redefine_cleaned = majiq_redefine_cleaned.reset_index(drop=True)
majiq_redefine_cleaned['exon_list'] = majiq_redefine_cleaned['exon_list'].apply(str)
majiq_redefine_cleaned['graph'] = majiq_redefine_cleaned['graph'].apply(str)

majiq_redefine_cleaned['exon_list'] = majiq_redefine_cleaned['exon_list'].apply(ast.literal_eval)
majiq_redefine_cleaned['graph'] = majiq_redefine_cleaned['graph'].apply(ast.literal_eval)

majiq_redefine_cleaned['split_id'] = majiq_redefine_cleaned['ID'] 
split_id = majiq_redefine_cleaned['split_id'].unique()
split_parts = np.array_split(split_id, 50)
majiq_part_index = 1
for sp in split_parts:
    sub_df = majiq_redefine_cleaned[majiq_redefine_cleaned['split_id'].isin(sp)]
    sub_df.to_csv( args.outdir + '/harm_majiq_p' + str(majiq_part_index) + '.csv')
    majiq_part_index += 1

print('majiq junction mapping done ....')
missing_gene = set(majiq['gene_name']) - set(majiq_redefine_cleaned['gene_name'])
print('the number of leafcutter total conatins {}'.format(len(set(majiq['gene_name']))))
print('the number of missing genes by leafcutter mapping junction library is {}'.format(len(missing_gene)))
print(missing_gene)