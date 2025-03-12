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
from itertools import chain

parser = argparse.ArgumentParser()
parser.add_argument('-infile', dest="infile",  help="csv path")
parser.add_argument('-outdir', dest="outdir", help="file path")
parser.add_argument('-Altbam', dest='Altbam', type=str, help="file path list")
parser.add_argument('-Refbam', dest='Refbam', type=str, help="file path list")
parser.add_argument('-inputjson', dest='inputjson', type=str, help="file path list")
parser.add_argument('-FCthreshold', dest='FCthreshold', type=float,  default=1.2, help='Fold change threshold (default = 1.2)')
parser.add_argument('-maxcount', dest='maxcount', type=float,  default=25, help='Fold change threshold (default = 25)')
args = parser.parse_args()

fc_thred, invfc_thred = args.FCthreshold, 1/args.FCthreshold

# print(fc_thred)
# print(invfc_thred)
df_name = args.infile
# print(args.Altbam)
# print(args.Refbam)

libdepths = load_json(args.inputjson)

try:
    alt_bam_files = ast.literal_eval(args.Altbam)
    ref_bam_files = ast.literal_eval(args.Refbam)

except ValueError:
    print("Invalid format. Please provide a valid list of file paths.")
    print(alt_bam_files)
except SyntaxError:
    print("Syntax error. Please provide a proper list format.")

df = pd.read_csv(df_name)
# df = junctioncounts_df_new(df, alt_bam_files, ref_bam_files, libdepths['alt_depth'], libdepths['ref_depth'])


try:
    if 'all_cleaned.csv' in df_name:
        outdir_folder = args.outdir.split("out")[0] + 'out/mxe/'
        df['annotation'] = 'mxe'
        df['splice'] = 'rmats'
        
    else:
        outdir_folder = args.outdir


    df = junctioncounts_df_countTolist(df, alt_bam_files, ref_bam_files)

    df['alt_count13'] = df['alt_count13'].apply(lambda x: pairwise_division(x, libdepths['alt_depth'].values()))
    df['alt_count23'] = df['alt_count23'].apply(lambda x: pairwise_division(x, libdepths['alt_depth'].values()))
    df['alt_count12'] = df['alt_count12'].apply(lambda x: pairwise_division(x, libdepths['alt_depth'].values()))
    df['ref_count13'] = df['ref_count13'].apply(lambda x: pairwise_division(x, libdepths['ref_depth'].values()))
    df['ref_count23'] = df['ref_count23'].apply(lambda x: pairwise_division(x, libdepths['ref_depth'].values()))
    df['ref_count12'] = df['ref_count12'].apply(lambda x: pairwise_division(x, libdepths['ref_depth'].values()))

    df['alt_count13_mean'] =  df['alt_count13'].apply(lambda x: sum(x) / len(x) if x else 0)
    df['alt_count12_mean'] =  df['alt_count12'].apply(lambda x: sum(x) / len(x) if x else 0)
    df['alt_count23_mean'] =  df['alt_count23'].apply(lambda x: sum(x) / len(x) if x else 0)
    df['ref_count13_mean'] =  df['ref_count13'].apply(lambda x: sum(x) / len(x) if x else 0)
    df['ref_count12_mean'] =  df['ref_count12'].apply(lambda x: sum(x) / len(x) if x else 0)
    df['ref_count23_mean'] =  df['ref_count23'].apply(lambda x: sum(x) / len(x) if x else 0)


    IRin_mask = (df['annotation'] == 'IR') & (df['splice'] == 'in') & (df['alt_count13_mean'] <= invfc_thred * df['ref_count13_mean']) & (df['ref_count13_mean'] > 0 )
    IRout_mask = (df['annotation'] == 'IR') & (df['splice'] == 'out') & (df['alt_count13_mean'] >= fc_thred * df['ref_count13_mean']) & (df['alt_count13_mean'] >0)


    ESin_mask = (
        (df['annotation'] == 'ES') & 
        (df['splice'] == 'in') & 
        (((df['alt_count13_mean'] <= invfc_thred * df['ref_count13_mean']) &  (df['ref_count13_mean'] > 0)).astype(int) + 
            ((df['alt_count12_mean'] >= fc_thred * df['ref_count12_mean']) & (df['alt_count12_mean'] >0)).astype(int) + 
            ((df['alt_count23_mean'] >= fc_thred * df['ref_count23_mean'])  & (df['alt_count23_mean'] >0)).astype(int) >= 2) 
    )
    ESout_mask = (
        (df['annotation'] == 'ES') & 
        (df['splice'] == 'out') & 
        (((df['alt_count13_mean'] >= fc_thred * df['ref_count13_mean'] ) & (df['alt_count13_mean'] >0) ).astype(int) + 
            ((df['alt_count12_mean'] <= invfc_thred * df['ref_count12_mean']) & (df['ref_count12_mean'] > 0)).astype(int) + 
            ((df['alt_count23_mean'] <= invfc_thred * df['ref_count23_mean']) & (df['ref_count23_mean'] >0)).astype(int) >= 2)
    )

    # ESin_mask =  (df['annotation'] == 'ES') & (df['splice'] == 'in') & (df['ESin_conditions_met'] >= 2)
    # ESout_mask =  (df['annotation'] == 'ES') & (df['splice'] == 'out')& (df['ESout_conditions_met'] >= 2)


    A3SSshortP_mask =  ((df['strand'] == '+' ) & (df['annotation'] == 'A3SS') 
                        &  (df['splice'] == 'short') 
                        & (df['alt_count13_mean'] >= fc_thred * df['ref_count13_mean']) 
                        & (df['alt_count12_mean'] <= invfc_thred * df['ref_count12_mean'])
                        & ((df['alt_count13_mean'] >0 ) & (df['ref_count12_mean'] >0 ))
                        )
    A5SSshortN_mask =  ((df['strand'] == '-' ) & (df['annotation'] == 'A5SS') 
                        & (df['splice'] == 'short') 
                        & (df['alt_count13_mean'] >= fc_thred * df['ref_count13_mean']) 
                        & (df['alt_count12_mean'] <= invfc_thred * df['ref_count12_mean']) 
                        & ((df['alt_count13_mean'] >0 ) & (df['ref_count12_mean'] >0 ))
                        )

    A3SSlongP_mask =  ((df['strand'] == '+' ) & (df['annotation'] == 'A3SS') 
                        & (df['splice'] == 'long') 
                        & (df['alt_count13_mean'] <= invfc_thred * df['ref_count13_mean']) 
                        & (df['alt_count12_mean'] >= fc_thred * df['ref_count12_mean'])
                        & ((df['ref_count13_mean'] >0 ) & (df['alt_count12_mean'] >0 ))
                        )

    A5SSlongN_mask =  ((df['strand'] == '-' ) & (df['annotation'] == 'A5SS') 
                        & (df['splice'] == 'long') 
                        & (df['alt_count13_mean'] <= invfc_thred * df['ref_count13_mean']) 
                        & (df['alt_count12_mean'] >= fc_thred* df['ref_count12_mean'])
                        & ((df['ref_count13_mean'] >0 ) & (df['alt_count12_mean'] >0 ))
                        )

    A5SSshortP_mask =  ((df['strand'] == '+' ) & (df['annotation'] == 'A5SS') 
                        & (df['splice'] == 'short') 
                        & (df['alt_count13_mean'] >= fc_thred * df['ref_count13_mean']) 
                        & (df['alt_count23_mean'] <= invfc_thred * df['ref_count23_mean']) 
                        & ((df['alt_count13_mean'] >0 ) & (df['ref_count23_mean'] >0 ))
                        )

    A3SSshortN_mask =  ((df['strand'] == '-' ) & (df['annotation'] == 'A3SS') 
                        & (df['splice'] == 'short') 
                        & (df['alt_count13_mean'] >= fc_thred * df['ref_count13_mean']) 
                        & (df['alt_count23_mean'] <= invfc_thred * df['ref_count23_mean']) 
                        & ((df['alt_count13_mean'] >0 ) & (df['ref_count23_mean'] >0 ))
                        )

    A5SSlongP_mask =  ((df['strand'] == '+' ) & (df['annotation'] == 'A5SS') 
                        & (df['splice'] == 'long') 
                        & (df['alt_count13_mean'] <= invfc_thred * df['ref_count13_mean']) 
                        & (df['alt_count23_mean'] >= fc_thred * df['ref_count23_mean']) 
                        & ((df['ref_count13_mean'] >0 ) & (df['alt_count23_mean'] >0 ))
                        )

    A3SSlongN_mask =  ((df['strand'] == '-' ) & (df['annotation'] == 'A3SS') 
                        & (df['splice'] == 'long') 
                        & (df['alt_count13_mean'] <= invfc_thred * df['ref_count13_mean']) 
                        & (df['alt_count23_mean'] >= fc_thred * df['ref_count23_mean']) 
                        & ((df['ref_count13_mean'] >0 ) & (df['alt_count23_mean'] >0 ))
                        )

    final_mask = (
        IRin_mask | IRout_mask | 
        ESin_mask | ESout_mask | 
        A3SSshortP_mask | A5SSshortN_mask | 
        A3SSlongP_mask | A5SSlongN_mask | 
        A5SSshortP_mask | A3SSshortN_mask | 
        A5SSlongP_mask | A3SSlongN_mask
    )
    df['junctionLabel'] = False  # Initialize all rows to 'False'
    df.loc[final_mask, 'junctionLabel'] = True


    df['ttest_pvalue'] = df.apply(lambda x: [perform_t_test(x['ref_count13'], x['alt_count13']), 
                                            perform_t_test(x['ref_count12'], x['alt_count12']),
                                            perform_t_test(x['ref_count23'], x['alt_count23'])], axis=1)

    # Apply the function to each list in the column 'ttest_pvalue'
    df['ttest_pvalue_removeNAs'] = df['ttest_pvalue'].apply(remove_nans)
    df['fisher_pvalue'] = df['ttest_pvalue_removeNAs'].apply(FisherP)
    try:
        #df['count_normalize'] = df[['alt_count13', 'alt_count12', 'alt_count23', 'ref_count13', 'ref_count12', 'ref_count23']].max(axis=1)
        df['count_normalize'] = df[['alt_count13', 'alt_count12', 'alt_count23', 'ref_count13', 'ref_count12', 'ref_count23']].apply(max_in_row, axis=1)
        val = args.maxcount/np.mean(list(libdepths['alt_depth'].values()) + list(libdepths['ref_depth'].values()))
        df['count_normalize_label'] = df['count_normalize'] >= val
    #df['count_normalize'] = df[['alt_count13', 'alt_count12', 'alt_count23', 'ref_count13', 'ref_count12', 'ref_count23']].apply(max_in_row, axis=1)
    except ValueError:
        print('count_normalize is NA: {}'.format(df_name))


    base_name = df_name.split('/')[-1][:-4]


    df.to_csv(outdir_folder + '/' + base_name + '_add_junction.csv', index=False)

except:
    print('Debug the script with: {}'.format(df_name))