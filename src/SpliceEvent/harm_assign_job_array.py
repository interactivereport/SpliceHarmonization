import numpy as np
import itertools
import pandas as pd
import warnings
import re
import glob
import argparse
warnings.filterwarnings('ignore')
from itertools import combinations
import os
from utils import *
#from mpi4py import MPI
import ast
# %%
def process_sj_file(csv_file, output_file, relabel_file, majiq_cutoff_val_file):
    print('is processing {}...'.format(csv_file))
    relabel = pd.read_csv(relabel_file)
    #introns_new = pd.read_csv(args.intron)
    sj_file = pd.read_csv(csv_file)
    sj_file['graph_str'] = sj_file['graph'].astype(str)
    clusters = sj_file.groupby(['ID', 'comparison', 'graph_str'])

    if 'harm_rmats_p' in csv_file:
        #clusters = sj_file['ID'].unique()
        file_clean = []
        #clusters = sj_file.groupby('ID')

        for _, sub_df in clusters:
            cassttee_groups =  [list(x) for x in set(tuple(ast.literal_eval(x)) for x in sub_df['casstte_list'])]
            relabel_dict = {}
            relabel_filtered = relabel[relabel['gene_name'] == sub_df['gene_name'].unique()[0]]
            for _, row in relabel_filtered.iterrows():
                exon_number = row['unified_exon']
                exon_coordinates = row[['Start', 'End']].to_list()
                relabel_dict[exon_number] = exon_coordinates
            for i, casstte in enumerate(cassttee_groups):
                l1to3r = [np.nan]*6
                for j, exon_ in enumerate(casstte):
                    l1to3r[j*2] =  relabel_dict[exon_][0]
                    l1to3r[j*2+1] =  relabel_dict[exon_][1]
                if len(casstte) < 3:
                    s_ = s_ = sub_df[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in casstte))]
                    s_ = s_[s_['exon_list'].apply(lambda x: str(x) == str(casstte))]
                    sub_score = float(s_['score'].iloc[0])
                    sub_PSI_1 = float(s_['PSI_1'].iloc[0])
                    sub_PSI_2 = float(s_['PSI_2'].iloc[0])
                    sub_exon_list = [s_['exon_list'].iloc[0]]
                else:
                    sub_score, sub_PSI_1, sub_PSI_2, sub_exon_list =  [], [], [], []
                        # print(sub_df['gene_name'].unique()[0])
                        # print()
                    for short_casstte in [casstte[:-1], casstte[1:]]:
                        #print(short_casstte)
                        s_ = sub_df[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte))]
                        s_ = s_[s_['casstte_list'].apply(lambda x: str(x) == str(casstte))]
                        if len(s_) >0:
                            sub_score_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'score']
                            sub_score.append([float(value) if value != 'NA' else np.nan for value in sub_score_series][0])
                            sub_PSI_1_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'PSI_1']
                            sub_PSI_1.append([float(value) if value != 'NA' else np.nan for value in sub_PSI_1_series][0])
                            sub_PSI_2_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'PSI_2']
                            sub_PSI_2.append([float(value) if value != 'NA' else np.nan for value in sub_PSI_2_series][0])
                            sub_exon_list.append(list(sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'exon_list'])[0])
                            
                sub_ID = sub_df.iloc[i]['ID']
                sub_comparison= sub_df.iloc[i]['comparison']
                sub_method= sub_df.iloc[i]['method']
                sub_gene_name= sub_df.iloc[i]['gene_name']
                sub_graph = sub_df.iloc[i]['graph']

                try:
                    #print(sub_score)
                    # print(list(sub_PSI_1))
                    # print(sub_score)
                    row_ = [sub_ID, sub_comparison, sub_method, sub_gene_name, sub_graph, np.nanmean(sub_score), np.nanmean(sub_PSI_1), np.nanmean(sub_PSI_2) , sub_exon_list,  casstte, relabel_filtered['Chromosome'].unique()[0], relabel_filtered['Strand'].unique()[0] ] + l1to3r
                    file_clean.append(row_) 
                except:
                    print('{} with graph{}'.fomart(sub_gene_name, sub_graph))

        file_clean = pd.DataFrame(file_clean, columns=['ID', 'comparison', 'method','gene_name',  'graph', 'score', 'PSI_1', 'PSI_2', 'exon_list',
                                            'cassette_list', 'chr', 'strand', 'cassette_1L', 'cassette_1R','cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R'])
        file_clean['cassette_list'] = file_clean['cassette_list'].astype(str)
        # file_clean = file_clean[file_clean.apply(lambda row: set(row['cassette_list']).issubset(row['graph']), axis=1)]
        file_clean = file_clean[file_clean.apply(lambda row: ((str(row['graph']) == str(row['cassette_list'])) and not ('MXE' in row['ID'])) or
                                                 (all(elem in ast.literal_eval(row['graph']) for elem in ast.literal_eval(row['cassette_list'])) and 'MXE' in row['ID']),
                                                 axis=1)]
        file_clean.to_csv(output_file)

    elif 'harm_leafcutter_p' in csv_file:
        file_clean = []
        for key_, sub_df in clusters:
            cassttee_groups =  [list(x) for x in set(tuple(ast.literal_eval(x)) for x in sub_df['casstte_list'])]
            relabel_dict = {}
            if len(sub_df['gene_name'].unique()) == 1:
                relabel_filtered = relabel[relabel['gene_name'] == sub_df['gene_name'].unique()[0]]
                for _, row in relabel_filtered.iterrows():
                    exon_number = row['unified_exon']
                    exon_coordinates = row[['Start', 'End']].to_list()
                    relabel_dict[exon_number] = exon_coordinates
                for i, casstte in enumerate(cassttee_groups):
                    l1to3r = [np.nan]*6
                    for j, exon_ in enumerate(casstte):
                        l1to3r[j*2] =  relabel_dict[exon_][0]
                        l1to3r[j*2+1] =  relabel_dict[exon_][1]
                    if len(casstte) < 3:
                        s_ = sub_df[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in casstte))]
                        s_ = s_[s_['exon_list'].apply(lambda x: str(x) == str(casstte))]
                        try:
                            sub_score = float(s_['score'].iloc[0])
                            sub_PSI_1 = float(s_['PSI_1'].iloc[0])
                            sub_PSI_2 = float(s_['PSI_2'].iloc[0])
                            sub_exon_list = [s_['exon_list'].iloc[0]]
                        except:
                            print('---line105 new debug:')
                            print(csv_file)
                            print(key_)
                            print(casstte)
                            print(s_)

                    else:
                        sub_score, sub_PSI_1, sub_PSI_2, sub_exon_list =  [], [], [], []
                        long_casstte = [casstte[0], casstte[-1]]
                        # s_long = sub_df[sub_df['exon_list'].apply(lambda x: all(str(elem) in x for elem in long_casstte))]
                        s_long = sub_df[sub_df['exon_list'] == str(long_casstte)]
                        if len(s_long) >0:
                            # print(sub_df['gene_name'].unique()[0])
                            # print()
                            for short_casstte in [casstte[:-1], casstte[1:]]:
                                #print(short_casstte)
                                #s_ = sub_df[sub_df['exon_list'].apply(lambda x: all(str(elem) in x for elem in short_casstte))]
                                s_ = sub_df[sub_df['exon_list'].apply(lambda x: str(x) == str(short_casstte))]
                                s_ = s_[s_['casstte_list'].apply(lambda x: str(x) == str(casstte))]
                                if len(s_) >0:
                                    sub_score_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'score']
                                    sub_score.append(max([float(value) if value != 'NA' else np.nan for value in sub_score_series][0], s_long['score'].to_list()[0]))
                                    sub_PSI_1_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'PSI_1']
                                    sub_PSI_1.append([float(value) if value != 'NA' else np.nan for value in sub_PSI_1_series][0]/ 
                                                    ([float(value) if value != 'NA' else np.nan for value in sub_PSI_1_series][0] + max(s_long['PSI_1'].to_list()[0], 1e-10)))
                                    sub_PSI_2_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'PSI_2']
                                    sub_PSI_2.append([float(value) if value != 'NA' else np.nan for value in sub_PSI_2_series][0]/
                                                    ([float(value) if value != 'NA' else np.nan for value in sub_PSI_2_series][0] + max(s_long['PSI_2'].to_list()[0], 1e-10)))
                                    sub_exon_list.append(list(sub_df.loc[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'exon_list'])[0])

                    sub_ID = sub_df.iloc[i]['ID']
                    sub_comparison= sub_df.iloc[i]['comparison']
                    sub_method= sub_df.iloc[i]['method']
                    sub_gene_name= sub_df.iloc[i]['gene_name']
                    sub_graph = sub_df.iloc[i]['graph']

                    try:
                        #print(sub_score)
                        # print(list(sub_PSI_1))
                        # print(sub_score)
                        row_ = [sub_ID, sub_comparison, sub_method, sub_gene_name, sub_graph, np.nanmean(sub_score), np.nanmean(sub_PSI_1), np.nanmean(sub_PSI_2) , sub_exon_list,  casstte, relabel_filtered['Chromosome'].unique()[0], relabel_filtered['Strand'].unique()[0] ] + l1to3r
                        file_clean.append(row_)
                    except:
                        print('{} with graph{}'.format(sub_gene_name, sub_graph))
            else:
                print('Cluster ID is{}, with Cluster gene name is {}'.format(sub_df['ID'].unique()[0], sub_df['gene_name'].unique()))
        file_clean = pd.DataFrame(file_clean, columns=['ID', 'comparison', 'method','gene_name',  'graph', 'score', 'PSI_1', 'PSI_2', 'exon_list',
                                            'cassette_list', 'chr', 'strand', 'cassette_1L', 'cassette_1R','cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R'])
        file_clean['cassette_list'] = file_clean['cassette_list'].astype(str)
        file_clean = file_clean[ ~ ((file_clean['PSI_1'] == 0.5) & (file_clean['PSI_2'] == 0.5))].reset_index(drop=True)
        file_clean.to_csv(output_file)

    else:
        majiq_cutoff_val = pd.read_csv(majiq_cutoff_val_file)
        file_clean = []
        for key_, sub_df in clusters:
            ###-----20250309 added
            sub_df = sub_df[['ID', 'score', 'comparison', 'method', 'PSI_1', 'PSI_2',
                    'type', 'graph', 'geneSymbol', 'start', 'end', 'chr', 'strand',
                    'coordinates', 'junction_filter', 'gene_name', 'junction_start',
                    'junction_end', 'unified_junction_label', 'exon_list', 'cutoff',
                    'casstte_list', 'repeated_exon', 'repeated_casstte', 'repeat_exon_len',
                    'repeat_casstee_len', 'split_id', 'graph_str',]].drop_duplicates().sort_values(by='score', ascending=False)
            sub_df = sub_df.drop_duplicates(subset= ['ID',  'comparison', 'method', 
                    'type', 'graph', 'geneSymbol', 'chr', 'strand', 'gene_name',  'unified_junction_label', 'exon_list', 'cutoff',
                    'casstte_list', 'repeated_exon', 'repeated_casstte', 'repeat_exon_len',
                    'repeat_casstee_len', 'split_id', 'graph_str'], keep='first')
            ###-----20250309 end
            cassttee_groups =  [list(x) for x in set(tuple(ast.literal_eval(x)) for x in sub_df['casstte_list'])]
            relabel_dict = {}
            sub_majiq = majiq_cutoff_val[majiq_cutoff_val['ID'] == key_[0]]
            ###-----20250309 added
            sub_majiq = sub_majiq[['score', 'comparison', 'method', 'PSI_1', 'PSI_2', 'ID', 'type',
            'graph', 'geneSymbol', 'start', 'end', 'chr', 'strand', 'coordinates',
            'junction_filter', 'gene_name', 'junction_start', 'junction_end',
            'unified_junction_label', 'exon_list', 'cutoff']].drop_duplicates().sort_values(by=['score', 'cutoff'], ascending=False)
            sub_majiq = sub_majiq.drop_duplicates(subset=[ 'comparison', 'method',  'ID', 'type',
            'graph', 'geneSymbol',  'chr', 'strand', 'gene_name', 'unified_junction_label', 'exon_list', 'cutoff'], keep ='first' )
            ###------20250309 end
            if len(sub_df['gene_name'].unique()) == 1:
                relabel_filtered = relabel[relabel['gene_name'] == sub_df['gene_name'].unique()[0]]
                for _, row in relabel_filtered.iterrows():
                    exon_number = row['unified_exon']
                    exon_coordinates = row[['Start', 'End']].to_list()
                    relabel_dict[exon_number] = exon_coordinates
                for i, casstte in enumerate(cassttee_groups):
                    l1to3r = [np.nan]*6
                    for j, exon_ in enumerate(casstte):
                        l1to3r[j*2] =  relabel_dict[exon_][0]
                        l1to3r[j*2+1] =  relabel_dict[exon_][1]
                    if len(casstte) < 3:
                        s_ = sub_df[sub_df['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in casstte))]
                        if len(s_) != 1:
                            print('---line180 debug:')
                            print(csv_file)
                            print(key_)
                            print(casstte)
                            print(s_)
                        # s_ = s_[s_['casstte_list'].apply(lambda x: str(x) == str(casstte))]
                        #sub_score_series = list(sub_majiq.loc[sub_majiq['exon_list'].apply(lambda x: all(str(elem) in x for elem in short_casstte)), 'score'])
                        #sub_score = float(s_['score'].iloc[0])
                        #sub_score = sub_majiq[(sub_majiq['exon_list'].apply(lambda x: all(str(elem) in x for elem in casstte)))]
                        sub_score = sub_majiq[(sub_majiq['exon_list'].apply(lambda x: str(x) == str(casstte)))]
                        #s_ = s_[s_['casstte_list'].apply(lambda x: all(str(elem) in x for elem in casstte))]
                        if len(sub_score) > 4:
                            print('ATTENTION: debug for less then 3 casette!')
                            print(csv_file)
                            print(key_)
                            print('---')
                            # raise ValueError("Error: Not find all cutoff val.")
                        sub_score_cutoff1 =  float(sub_score['score'].iloc[0])
                        sub_score_cutoff2 =  float(sub_score['score'].iloc[1])
                        sub_score_cutoff3 =  float(sub_score['score'].iloc[2])
                        sub_score_cutoff4 =  float(sub_score['score'].iloc[3])
                        try:
                            sub_PSI_1 = float(s_['PSI_1'].iloc[0])
                            sub_PSI_2 = float(s_['PSI_2'].iloc[0])
                            sub_exon_list = [s_['exon_list'].iloc[0]]
                        except:
                            print('---line198 debug:')
                            print(csv_file)
                            print(key_)
                            print(casstte)
                            print(s_)
                    else:
                        sub_score_cutoff1, sub_score_cutoff2, sub_score_cutoff3, sub_score_cutoff4,  sub_PSI_1, sub_PSI_2, sub_exon_list = [],[],[], [], [], [], []
                        long_casstte = [casstte[0], casstte[-1]]
                        # s_long = sub_df[sub_df['exon_list'].apply(lambda x: all(str(elem) in x for elem in long_casstte))]
                        s_long = sub_df[sub_df['exon_list'] == str(long_casstte)]
                        s_long = s_long[s_long['casstte_list'].apply(lambda x: str(x) == str(casstte))]
                        s_long_cutoff_score = sub_majiq.loc[sub_majiq['exon_list'] == str(long_casstte) , 'score']
                        sub_long_cutoff_all = [float(value) if value != 'NA' else np.nan for value in s_long_cutoff_score]
                        if len(s_long) == 1:
                            # print(sub_df['gene_name'].unique()[0])
                            # print()
                            for short_casstte in [casstte[:-1], casstte[1:]]:
                                #print(short_casstte)
                                #s_ = sub_df[sub_df['exon_list'].apply(lambda x: all(str(elem) in x for elem in short_casstte))]
                                s_ = sub_df[sub_df['exon_list'].apply(lambda x: str(x) == str(short_casstte))]
                                s_ = s_[s_['casstte_list'].apply(lambda x: str(x) == str(casstte))]
                                if len(s_) == 1:
                                    #sub_score_series = sub_df.loc[sub_df['exon_list'].apply(lambda x: all(str(elem) in x for elem in short_casstte)), 'score']
                                    sub_score_series = list(sub_majiq.loc[sub_majiq['exon_list'].apply(lambda x: str(x) == str(short_casstte )), 'score'])
                                    if len(sub_score_series) > 4:
                                        print(short_casstte)
                                        print(casstte)
                                        print("3 cassette fail")
                                        print(csv_file)
                                        print(key_)
                                        print('---')
                                        raise ValueError("Error: Not find all cutoff val.")
                                    
                                    sub_score_cutoff_all = [float(value) if value != 'NA' else np.nan for value in sub_score_series]
                                    

                                    #sub_score.append(max([float(value) if value != 'NA' else np.nan for value in sub_score_series][0], s_long['score'].to_list()[0]))
                                    try:
                                        sub_score_cutoff1.append(max(sub_score_cutoff_all[0], sub_long_cutoff_all[0]))
                                        sub_score_cutoff2.append(max(sub_score_cutoff_all[1], sub_long_cutoff_all[1]))
                                        sub_score_cutoff3.append(max(sub_score_cutoff_all[2], sub_long_cutoff_all[2]))
                                        sub_score_cutoff4.append(max(sub_score_cutoff_all[3], sub_long_cutoff_all[3]))
                                    except:
                                        print('---debug line 245')
                                        print(short_casstte)
                                        print(casstte)
                                        print(csv_file)
                                        print(key_)
                                        print('---')


                                    sub_PSI_1_series = s_.loc[s_['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'PSI_1']
                                    sub_PSI_1.append([float(value) if value != 'NA' else np.nan for value in sub_PSI_1_series][0]/ 
                                                    ([float(value) if value != 'NA' else np.nan for value in sub_PSI_1_series][0] + max(s_long['PSI_1'].to_list()[0], 1e-10)))
                                    sub_PSI_2_series = s_.loc[s_['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x)for elem in short_casstte)), 'PSI_2']
                                    sub_PSI_2.append([float(value) if value != 'NA' else np.nan for value in sub_PSI_2_series][0]/
                                                    ([float(value) if value != 'NA' else np.nan for value in sub_PSI_2_series][0] + max(s_long['PSI_2'].to_list()[0], 1e-10)))
                                    sub_exon_list.append(list(s_.loc[s_['exon_list'].apply(lambda x: all(elem in ast.literal_eval(x) for elem in short_casstte)), 'exon_list'])[0])
                                
                                elif len(s_) > 1:
                                    print(len(s_))
                                    print('ATTENTION: debug for 3 casette: short junction')
                                    print(short_casstte)
                                    print(key_)
                                    print(casstte)
                                    print('---')
                                else:
                                    continue
                        elif len(s_long) > 1:
                            print('ATTENTION: debug for 3 casette: long junction')
                            print(s_long)
                            print(key_)
                            print('---')
                        else:
                            continue

                    sub_ID = key_[0]
                    sub_comparison= sub_df.iloc[i]['comparison']
                    sub_method= sub_df.iloc[i]['method']
                    sub_gene_name= sub_df.iloc[i]['gene_name']
                    sub_graph = sub_df.iloc[i]['graph']

                    try:
                        #print(sub_score)
                        # print(list(sub_PSI_1))
                        # print(sub_score)
                        #row_ = [sub_ID, sub_comparison, sub_method, sub_gene_name, sub_graph, np.nanmean(sub_score), np.nanmean(sub_PSI_1), np.nanmean(sub_PSI_2) , sub_exon_list,  casstte, relabel_filtered['Chromosome'].unique()[0], relabel_filtered['Strand'].unique()[0] ] + l1to3r
                        row_ = [sub_ID, sub_comparison, sub_method, sub_gene_name, sub_graph, np.nanmean(sub_score_cutoff1), np.nanmean(sub_score_cutoff2), np.nanmean(sub_score_cutoff3), np.nanmean(sub_score_cutoff4),
                                np.nanmean(sub_PSI_1), np.nanmean(sub_PSI_2) , sub_exon_list,  casstte, relabel_filtered['Chromosome'].unique()[0], relabel_filtered['Strand'].unique()[0] ] + l1to3r
                        
                        file_clean.append(row_)
                    except:
                        print('{} with graph{}'.format(sub_gene_name, sub_graph))
            else:
                print('Cluster ID is{}, with Cluster gene name is {}'.format(sub_df['ID'].unique()[0], sub_df['gene_name'].unique()))
        file_clean = pd.DataFrame(file_clean, columns=['ID', 'comparison', 'method','gene_name',  'graph', 'score_cutoff0.1', 'score_cutoff0.2', 'score_cutoff0.3', 'score_cutoff0.4', 'PSI_1', 'PSI_2', 'exon_list',
                                            'cassette_list', 'chr', 'strand', 'cassette_1L', 'cassette_1R','cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R'])
        file_clean['cassette_list'] = file_clean['cassette_list'].astype(str)
        file_clean = file_clean[ ~((file_clean['PSI_1'] == 0.5) & (file_clean['PSI_2'] == 0.5) & ~(file_clean['ID'].str.contains('ir'))) ].reset_index(drop=True)
        file_clean.to_csv(output_file)

    print('{} procssing is done'.format(csv_file))
    missing_gene = set(sj_file['gene_name']) - set(file_clean['gene_name'])
    print('the number of total genes before mapping is {}'.format(len(set(sj_file['gene_name']))))
    print('the number of missing genes by mapping junction library is {}'.format(len(missing_gene)))
    print(missing_gene)
    print('---')
# %%
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-exonlabel', dest="exonlabel", help ='csv path')
    parser.add_argument('-infile', dest="infile",  help=".csv path")
    parser.add_argument('-outdir', dest="outdir", help="save  .csv_path")
    parser.add_argument('-majiqcutofffile', dest='majiqcutofffile', help= 'majiq cutoff file path')
    args = parser.parse_args()

    filename = args.infile
    relabel_file = args.exonlabel
    if 'mxe' in filename:
        output_file = os.path.join(os.path.dirname(args.outdir), 'mxe', 'out',  os.path.basename(filename).rstrip('.csv') + '_sj.output.csv' )
    else:
        output_file = os.path.join(args.outdir, os.path.basename(filename).rstrip('.csv') + '_sj.output.csv')
    # output_file = os.path.join(args.outdir, os.path.basename(filename).rstrip('.csv') + '_sj.output.csv')
    majiq_cutoff_val_file = args.majiqcutofffile
    process_sj_file(filename, output_file, relabel_file, majiq_cutoff_val_file)