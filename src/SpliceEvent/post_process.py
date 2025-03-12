import numpy as np
import itertools
import pandas as pd
import warnings
import re
warnings.filterwarnings('ignore')
from itertools import combinations
import ast
import argparse
import glob
import os
from utils import *
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('-indir', dest="indir",  help="csv path")
parser.add_argument('-outdir', dest="outdir", help="file path")
parser.add_argument('-allgene', dest="allgene", help="file path")
parser.add_argument('-majiq_score_cutoff', dest="majiq_score_cutoff",  type=float, default=0.9, help="majiq score cutoff for filter")


# parser.add_argument('-Altbam', dest='Altbam', type=str, help="file path list")
# parser.add_argument('-Refbam', dest='Refbam', type=str, help="file path list")
args = parser.parse_args()

WD = Path(args.outdir).resolve().parents[0]
majiq_junction_4IR = pd.read_csv(WD / 'junction_prep/majiq_junction_prep.csv')
majiq_junction_4IR = majiq_junction_4IR.loc[(majiq_junction_4IR['lsv_type']== 'i'), ]

rmats_files = glob.glob(args.indir + '/*rmats_*sj.output.csv')
leafcutter_files = glob.glob(args.indir + '/*leafcutter_*sj.output.csv')
majiq_files = glob.glob(args.indir + '/*majiq_*sj.output.csv')

# try:
#     alt_bam_files = ast.literal_eval(args.Altbam)
#     ref_bam_files = ast.literal_eval(args.Refbam)
#     print(alt_bam_files)
#     print(ref_bam_files)
# except ValueError:
#     print("Invalid format. Please provide a valid list of file paths.")
# except SyntaxError:
#     print("Syntax error. Please provide a proper list format.")

rmats = []
for f in rmats_files:
    df = pd.read_csv(f)
    df = df.drop(columns=['Unnamed: 0'])
    rmats.append(df)

leafcutters = []
for f in leafcutter_files:
    df = pd.read_csv(f)
    df = df.drop(columns=['Unnamed: 0'])
    leafcutters.append(df) 

majiqs = []
for f in majiq_files:
    df = pd.read_csv(f)
    df = df.drop(columns=['Unnamed: 0'])
    majiqs.append(df)

majiq = pd.concat(majiqs)
try:
    rmat= pd.concat(rmats)
    rmat['dpsi'] = rmat['PSI_2'] - rmat['PSI_1']

    rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_3L'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2L']
    rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_3R'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2R']
    rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2L'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_1L']
    no_rmat_flag = 0
except:
    print('no rmats result')
    no_rmat_flag = 1

try: 
    leafcutter = pd.concat(leafcutters)
    leafcutter['dpsi'] = leafcutter['PSI_2'] - leafcutter['PSI_1']
    leafcutter_IR_candidate = leafcutter[leafcutter['cassette_3L'].isna()].reset_index(drop=True)
    leafcutter_non_IR = leafcutter[~leafcutter['cassette_3L'].isna()].reset_index(drop=True)
    leafcutter_IR_candidate['graph_len'] = leafcutter_IR_candidate['graph'].apply(lambda x: len(ast.literal_eval(x)))
    leafcutter_IR_candidate.loc[(leafcutter_IR_candidate['graph_len'] == 2) & (leafcutter_IR_candidate['ID'].map(leafcutter_IR_candidate['ID'].value_counts() == 1)), 'cassette_3L'] = \
        leafcutter_IR_candidate.loc[(leafcutter_IR_candidate['graph_len'] == 2) & (leafcutter_IR_candidate['ID'].map(leafcutter_IR_candidate['ID'].value_counts() == 1)), 'cassette_2L']
    leafcutter_IR_candidate.loc[(leafcutter_IR_candidate['graph_len'] == 2) & (leafcutter_IR_candidate['ID'].map(leafcutter_IR_candidate['ID'].value_counts() == 1)), 'cassette_3R'] = \
        leafcutter_IR_candidate.loc[(leafcutter_IR_candidate['graph_len'] == 2) & (leafcutter_IR_candidate['ID'].map(leafcutter_IR_candidate['ID'].value_counts() == 1)), 'cassette_2R']
    leafcutter_IR_candidate.loc[(leafcutter_IR_candidate['graph_len'] == 2) & (leafcutter_IR_candidate['ID'].map(leafcutter_IR_candidate['ID'].value_counts() == 1)), 'cassette_2L'] = \
        leafcutter_IR_candidate.loc[(leafcutter_IR_candidate['graph_len'] == 2) & (leafcutter_IR_candidate['ID'].map(leafcutter_IR_candidate['ID'].value_counts() == 1)), 'cassette_1L']
    leafcutter_new = pd.concat([leafcutter_IR_candidate, leafcutter_non_IR])
    no_leafcutter_flag = 0 
except:
    print('no leafcutter result')
    no_leafcutter_flag = 1

# majiq['PdPSI'] = majiq['score']
# rmat['FDR'] = 1 - rmat['score']
# leafcutter['FDR'] = 1 - leafcutter['score']

# Apply the function to each row of the DataFrame and create new columns
majiq[['score_list', 'score', 'score_index']] = majiq.apply(lambda x: majiq_score_list(x, args.majiq_score_cutoff), axis=1)
majiq['dpsi'] = majiq['PSI_2'] - majiq['PSI_1']


# rmat['dpsi'] = rmat['PSI_2'] - rmat['PSI_1']

# rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_3L'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2L']
# rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_3R'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2R']
# rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2L'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_1L']

#################### 
# DO IT LATER IN JUNCTION FILTER PART
# rmat_MXE = rmat[rmat['ID'].str.contains('MXE')]
# rmat_MXE_candidate = rmat.loc[ (rmat['ID'].str.contains('MXE')) ]
# rmat_non_MXE = rmat.loc[~rmat.index.isin(rmat_MXE_candidate.index)]

# rmat_MXE_candidate[['e1', 'e2', 'e3']] = pd.DataFrame(
#     rmat_MXE_candidate['cassette_list'].apply(ast.literal_eval).tolist(),
#     index=rmat_MXE_candidate.index)
# rmat_MXE_candidate['graph'] = rmat_MXE_candidate['graph'].apply(ast.literal_eval)
# rmat_MXE_candidate = rmat_MXE_candidate.loc[rmat_MXE_candidate.apply(lambda x: (x['e1'] == x['graph'][0]) & (x['e3'] == x['graph'][-1]), axis=1)]
# rmat_MXE_candidate['graph'] = rmat_MXE_candidate['graph'].astype(str)
# rmat_MXE_candidate = rmat_MXE_candidate.groupby(['gene_name', 'graph', 'e1', 'e3']).filter(lambda x: x['e2'].nunique() == 2)
# rmat_MXE_candidate['MXEID']  =  rmat_MXE_candidate['ID'] + ':' + rmat_MXE_candidate['gene_name'] + ':'  + rmat_MXE_candidate['graph']
# rmat_MXE_candidate = rmat_MXE_candidate.drop(columns=['e1', 'e2', 'e3'])
# rmat_new = pd.concat([rmat_MXE_candidate, rmat_non_MXE])
####################

# majiq['cutoff'] = majiq['comparison'].apply(lambda x: np.float64(x.split('_')[1]))
# sig = majiq.groupby(['gene_name', 'cassette_list']).filter(lambda x: x['score'].max() >= 0.95)
# sig = sig.reset_index(drop=True)
# nosig = majiq.groupby(['gene_name', 'cassette_list']).filter(lambda x: x['score'].max() < 0.95)
# nosig = nosig.reset_index(drop=True)

# sig_new = sig.loc[sig.reset_index().groupby(['gene_name', 'cassette_list'])['cutoff'].idxmax()].reset_index(drop=True)
# nosig_new = nosig.loc[nosig.reset_index().groupby(['gene_name', 'cassette_list'])['cutoff'].idxmin()].reset_index(drop=True)

# overlap = pd.merge(sig_new, nosig_new, on=['gene_name', 'cassette_list'], how='inner')

# # If there's overlap, remove those rows from df_2
# if not overlap.empty:
#     # Using isin to filter out overlapping rows in df_2
#     condition = (sig_new['gene_name'].isin(overlap['gene_name'])) & (sig_new['cassette_list'].isin(overlap['cassette_list']))
#     nosig_new = nosig_new[~condition]

# majiq = pd.concat([sig_new, nosig_new]).drop(columns=['cutoff'])


majiq_IR_candidate = majiq.loc[ (majiq['ID'].str.contains('ir')) & (majiq['cassette_list'].apply(lambda x: len(ast.literal_eval(x)) == 2))]
majiq_non_IR = majiq.loc[~majiq.index.isin(majiq_IR_candidate.index)]
majiq_IR_candidate['cassette_3L'] = majiq_IR_candidate['cassette_2L']
majiq_IR_candidate['cassette_3R'] = majiq_IR_candidate['cassette_2R']
majiq_IR_candidate['cassette_2L'] = majiq_IR_candidate['cassette_1L']
majiq_IR_candidate =  pd.merge(majiq_IR_candidate, majiq_junction_4IR[['start', 'end']], how='inner', left_on='cassette_1R', right_on='start')
majiq_IR_candidate = majiq_IR_candidate[(majiq_IR_candidate['cassette_3L'] + 1 == majiq_IR_candidate['end'] )]
majiq_IR_candidate['majiq_ir_annotation'] = True
majiq_IR_candidate.drop(columns=['start', 'end'])
majiq_new = pd.concat([majiq_IR_candidate, majiq_non_IR])
# majiq_IR_candidate = majiq[majiq['cassette_3L'].isna()].reset_index(drop=True)
# majiq_non_IR = majiq[~majiq['cassette_3L'].isna()].reset_index(drop=True)
# majiq_IR_candidate['graph_len'] = majiq_IR_candidate['graph'].apply(lambda x: len(ast.literal_eval(x)))
# majiq_IR_candidate.loc[(majiq_IR_candidate['graph_len'] == 2) & (majiq_IR_candidate['ID'].map(majiq_IR_candidate['ID'].value_counts() == 1)), 'cassette_3L'] = \
#     majiq_IR_candidate.loc[(majiq_IR_candidate['graph_len'] == 2) & (majiq_IR_candidate['ID'].map(majiq_IR_candidate['ID'].value_counts() == 1)), 'cassette_2L']
# majiq_IR_candidate.loc[(majiq_IR_candidate['graph_len'] == 2) & (majiq_IR_candidate['ID'].map(majiq_IR_candidate['ID'].value_counts() == 1)), 'cassette_3R'] = \
#     majiq_IR_candidate.loc[(majiq_IR_candidate['graph_len'] == 2) & (majiq_IR_candidate['ID'].map(majiq_IR_candidate['ID'].value_counts() == 1)), 'cassette_2R']
# majiq_IR_candidate.loc[(majiq_IR_candidate['graph_len'] == 2) & (majiq_IR_candidate['ID'].map(majiq_IR_candidate['ID'].value_counts() == 1)), 'cassette_2L'] = \
#     majiq_IR_candidate.loc[(majiq_IR_candidate['graph_len'] == 2) & (majiq_IR_candidate['ID'].map(majiq_IR_candidate['ID'].value_counts() == 1)), 'cassette_1L']
# majiq_new = pd.concat([majiq_IR_candidate, majiq_non_IR])

if (no_rmat_flag == 1) and (no_leafcutter_flag == 0):
    df = pd.concat([leafcutter_new, majiq_new])
elif (no_rmat_flag == 0) and (no_leafcutter_flag == 1):
    df = pd.concat([rmat, majiq_new])
elif (no_rmat_flag == 1) and (no_leafcutter_flag == 1):
    df = majiq_new
else:
    df = pd.concat([rmat, leafcutter_new, majiq_new])

df['eventID'] = df.groupby(['gene_name', 'cassette_list']).ngroup()
df['eventID'] = df['gene_name'] + '_' + df['eventID'].astype(str)

df[['cassette_3L', 'cassette_3R']] = df[['cassette_3L', 'cassette_3R']].astype('Int64')

df_pre = df 
df = df.dropna(subset=['cassette_1L', 'cassette_1R',
                   'cassette_2L', 'cassette_2R', 
                   'cassette_3L', 'cassette_3R' ])

# df = junctioncounts_df_new(df, alt_bam_files, ref_bam_files)

df.to_csv(args.outdir + '/all_cleaned.csv', index=False)

missing_gene = (set(df_pre['gene_name']) - set(df['gene_name']))
print('the number of total conatins {}'.format(len(set(df_pre['gene_name']))))
print('the number of missing genes is {}'.format(len(missing_gene)))
print(missing_gene)

# df  = df.dropna(subset=['cassette_3L'])
genes = df['gene_name'].unique()
for g_ in genes:
    geneDf = df[df['gene_name'] == g_]
    geneDf = geneDf.reset_index(drop=True)
    mask_ES = (geneDf['cassette_2L'] > geneDf['cassette_1R']) & (geneDf['cassette_2R'] < geneDf['cassette_3L']) & (geneDf['method'] != 'rmats')
    mask_a5ss_p = (geneDf['cassette_2L'] == geneDf['cassette_1L']) & (geneDf['cassette_2R'] != geneDf['cassette_3R']) & (geneDf['strand'] == '+') & (geneDf['method'] != 'rmats')
    mask_a3ss_n = (geneDf['cassette_2L'] == geneDf['cassette_1L']) & (geneDf['cassette_2R'] != geneDf['cassette_3R']) & (geneDf['strand'] == '-') & (geneDf['method'] != 'rmats')
    mask_a3ss_p = (geneDf['cassette_2R'] == geneDf['cassette_3R']) & (geneDf['cassette_2L'] != geneDf['cassette_1L']) & (geneDf['strand'] == '+') & (geneDf['method'] != 'rmats')
    mask_a5ss_n = (geneDf['cassette_2R'] == geneDf['cassette_3R']) & (geneDf['cassette_2L'] != geneDf['cassette_1L']) & (geneDf['strand'] == '-') & (geneDf['method'] != 'rmats')
    mask_drop_L = (geneDf['cassette_2R'] == geneDf['cassette_1R']) & (geneDf['method'] != 'rmats')
    mask_drop_R = (geneDf['cassette_3L'] == geneDf['cassette_2L']) & (geneDf['method'] != 'rmats')
    mask_IR = (geneDf['cassette_2R'] == geneDf['cassette_3R']) & (geneDf['cassette_2L'] == geneDf['cassette_1L']) 
    geneDf.loc[geneDf['method'] == 'rmats', 'annotation_rmats'] = geneDf.loc[geneDf['method'] == 'rmats', 'ID'].apply(lambda row: 
                                                                                                                'ES' if 'ES' in row 
                                                                                                                else 'MXE' if 'MXE' in row 
                                                                                                                else 'IR' if 'IR' in row 
                                                                                                                else 'A5SS' if 'A5SS' in row 
                                                                                                                else 'A3SS' if 'A3SS' in row
                                                                                                                else '')

    mask_ES = (geneDf['cassette_2L'] > geneDf['cassette_1R']) & (geneDf['cassette_2R'] < geneDf['cassette_3L']) 
    mask_a5ss_p = (geneDf['cassette_2L'] == geneDf['cassette_1L']) & (geneDf['cassette_2R'] != geneDf['cassette_3R']) & (geneDf['strand'] == '+') 
    mask_a3ss_n = (geneDf['cassette_2L'] == geneDf['cassette_1L']) & (geneDf['cassette_2R'] != geneDf['cassette_3R']) & (geneDf['strand'] == '-') 
    mask_a3ss_p = (geneDf['cassette_2R'] == geneDf['cassette_3R']) & (geneDf['cassette_2L'] != geneDf['cassette_1L']) & (geneDf['strand'] == '+') 
    mask_a5ss_n = (geneDf['cassette_2R'] == geneDf['cassette_3R']) & (geneDf['cassette_2L'] != geneDf['cassette_1L']) & (geneDf['strand'] == '-') 
    mask_drop_L = (geneDf['cassette_2R'] == geneDf['cassette_1R']) & (geneDf['method'] != 'rmats')
    mask_drop_R = (geneDf['cassette_3L'] == geneDf['cassette_2L']) & (geneDf['method'] != 'rmats')
    mask_IR = (geneDf['cassette_2R'] == geneDf['cassette_3R']) & (geneDf['cassette_2L'] == geneDf['cassette_1L']) 

    geneDf.loc[mask_ES, 'annotation'] = 'ES'
    geneDf.loc[mask_a5ss_p, 'annotation'] = 'A5SS'
    geneDf.loc[mask_a3ss_p, 'annotation'] = 'A3SS'
    geneDf.loc[mask_a5ss_n, 'annotation'] = 'A5SS'
    geneDf.loc[mask_a3ss_n, 'annotation'] = 'A3SS'
    geneDf.loc[mask_IR, 'annotation'] = 'IR'
    #geneDf.loc[(geneDf['majiq_ir_annotation'] ==True), 'annotation'] = 'IR'

    geneDf = geneDf[~(geneDf['method'] == 'rmats') | ((geneDf['method'] == 'rmats') & (geneDf['annotation_rmats'] != 'MXE') & (geneDf['annotation'] == geneDf['annotation_rmats']))|
                    ((geneDf['method'] == 'rmats') & (geneDf['annotation_rmats'] == 'MXE'))]

    try:

        geneDf.loc[(((geneDf['annotation'] == 'ES' ) | (geneDf['annotation'] == 'MXE' ) | (geneDf['annotation'] == 'IR' ) ) & (geneDf['dpsi'] >0 )), 'splice'] = 'in'
        geneDf.loc[(((geneDf['annotation'] == 'ES' ) | (geneDf['annotation'] == 'MXE' ) | (geneDf['annotation'] == 'IR' ) ) & (geneDf['dpsi'] <=0 )), 'splice'] = 'out'
        geneDf.loc[(((geneDf['annotation'] == 'A5SS') | (geneDf['annotation'] == 'A3SS')) & (geneDf['dpsi'] >0 )), 'splice'] = 'long'
        geneDf.loc[(((geneDf['annotation'] == 'A5SS') | (geneDf['annotation'] == 'A3SS')) & (geneDf['dpsi'] <= 0 )), 'splice'] = 'short'
        geneDf.drop(geneDf[mask_drop_L | mask_drop_R].index, inplace =True)
    
        # geneDf = geneDf[geneDf['annotation'] != 'MXE']
        sorted_genedf = geneDf.sort_values(by=['score', 'dpsi'], ascending=[False, False], key=lambda x: abs(x))
        sorted_genedf.reset_index(drop=True)
        sorted_genedf.to_csv(args.allgene + '/' + g_ + '.csv')
    
    except:
        print('...pls check ...')
        print(g_)