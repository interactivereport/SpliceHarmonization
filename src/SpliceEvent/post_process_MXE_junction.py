import pandas as pd
import pyranges as pr
import seaborn as sns
import sys
from utils import *
import matplotlib.pyplot as plt
import ast
import warnings
import argparse
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('-indir', dest="indir",  help="csv path")
parser.add_argument('-mxeindir', dest='mxeindir', help='csv path')
parser.add_argument('-outdir', dest="outdir", help="file path")
args = parser.parse_args()

df_current = []
for gene_file in glob.glob(args.indir + '/*.csv'):
    df_ = pd.read_csv(gene_file)
    df_current.append(df_)
df_current_concat = pd.concat(df_current)
df_current_concat['event_label'] =(df_current_concat['gene_name'] + ':' +
                        df_current_concat['cassette_1R'].astype(str) + ':' +
                        df_current_concat['cassette_2L'].astype(str) + ':' +
                        df_current_concat['cassette_2R'].astype(str) + ':' +
                        df_current_concat['cassette_3L'].astype(str))


mxe_rmats = pd.read_csv(args.mxeindir+'/all_cleaned_add_junction.csv')

rmats = pd.concat([df_current_concat[df_current_concat['method'] == 'rmats'], mxe_rmats])
non_rmats = df_current_concat[df_current_concat['method'] != 'rmats']

non_rmats[['e1', 'e2', 'e3']] = pd.DataFrame(
    non_rmats['cassette_list'].apply(ast.literal_eval).tolist(),
    index=non_rmats.index
)

filtered_non_rmats = non_rmats.groupby(['gene_name', 'ID', 'e1', 'e3']).filter(
    lambda x: (x['e2'].nunique() == 2) and (len(x['dpsi']) == 2) and (np.sign(x['dpsi'].iloc[0]) * np.sign(x['dpsi'].iloc[1]) < 0)
)

filtered_non_rmats['MXE_pseduoLabel'] = filtered_non_rmats['gene_name'] + ":"  + filtered_non_rmats['ID']  + ":" + filtered_non_rmats['e1'].astype(str) + ":" + filtered_non_rmats['e3'].astype(str)
filtered_non_rmats = filtered_non_rmats.sort_values('MXE_pseduoLabel')
non_rmats_result = filtered_non_rmats.groupby(['ID', 'gene_name', 'MXE_pseduoLabel', 'method', 'comparison', 'graph', 'chr', 'strand']).agg({
    "score": list,
    'exon_list':list,
    'cassette_list':unique_flattened_list,
    "cassette_1L": 'first',
    "cassette_1R": 'first',
    "cassette_2L": list,
    "cassette_2R": list,
    "cassette_3L": 'first',
    "cassette_3R": 'first',
    "junctionLabel": 'min',
    "count_normalize_label": 'min',
    'dpsi': list,
}).reset_index()

non_rmats_result = non_rmats_result.rename(columns={
    'score': 'score_list',
    'dpsi': 'dpsi_list'
})

non_rmats_result['score'] = non_rmats_result['score_list'].apply(lambda x: np.mean(x))
non_rmats_result['dpsi'] = non_rmats_result['dpsi_list'].apply(lambda x: np.mean(np.abs(x)))
non_rmats_result['splice'] = non_rmats_result['dpsi_list'].apply( lambda x: 'LinRout' if (x[0] > 0 and x[1] < 0) else 'RinLout' if (x[0] < 0 and x[1] > 0) else 'NA')
non_rmats_result['annotation'] = 'MXE'

rmats[['e1', 'e2', 'e3']] = pd.DataFrame(
    rmats['cassette_list'].apply(ast.literal_eval).tolist(),
    index=rmats.index
)

filtered_rmats = rmats.groupby(['gene_name', 'ID', 'e1', 'e3']).filter(
    lambda x: (x['e2'].nunique() == 2) 
)
filtered_rmats['MXE_pseduoLabel'] = filtered_rmats['gene_name'] + ":"  + filtered_rmats['ID']  + ":" + filtered_rmats['e1'].astype(str) + ":" + filtered_rmats['e3'].astype(str)
filtered_rmats = filtered_rmats.sort_values('MXE_pseduoLabel')
rmats_result = filtered_rmats.groupby(['ID', 'gene_name', 'MXE_pseduoLabel', 'method', 'comparison', 'graph', 'chr', 'strand']).agg({
    "score": list,
    'exon_list':list,
    'cassette_list':unique_flattened_list,
    "cassette_1L": 'first',
    "cassette_1R": 'first',
    "cassette_2L": list,
    "cassette_2R": list,
    "cassette_3L": 'first',
    "cassette_3R": 'first',
    "junctionLabel": 'min',
    "count_normalize_label": 'min',
    'dpsi': list,
}).reset_index()


rmats_result = rmats_result.rename(columns={
    'score': 'score_list',
    'dpsi': 'dpsi_list'
})

rmats_result['score'] = rmats_result['score_list'].apply(lambda x: np.mean(x))
rmats_result['dpsi'] = rmats_result['dpsi_list'].apply(lambda x: np.mean(x))
#rmats_result['splice'] = rmats_result['dpsi_list'].apply( lambda x: 'LtoR' if (x[0] > 0 and x[1] < 0) else 'RtoL' if (x[0] < 0 and x[1] > 0) else 'NA')
rmats_result['annotation'] = 'MXE'

mxe = pd.concat([non_rmats_result, rmats_result])
mxe.to_csv(args.outdir + '/MXE_result_add_junction.csv')
