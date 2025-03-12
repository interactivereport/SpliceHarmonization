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

# parser.add_argument('-Altbam', dest='Altbam', type=str, help="file path list")
# parser.add_argument('-Refbam', dest='Refbam', type=str, help="file path list")
args = parser.parse_args()

WD = Path(args.outdir).resolve().parents[0]

rmats_files = glob.glob(args.indir + '/*rmats_*sj.output.csv')

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

rmat= pd.concat(rmats)
# majiq['PdPSI'] = majiq['score']
# rmat['FDR'] = 1 - rmat['score']
# leafcutter['FDR'] = 1 - leafcutter['score']

rmat['dpsi'] = rmat['PSI_2'] - rmat['PSI_1']

rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_3L'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2L']
rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_3R'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2R']
rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_2L'] = rmat.loc[rmat['ID'].str.contains('IR'), 'cassette_1L']

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



df = rmat

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
