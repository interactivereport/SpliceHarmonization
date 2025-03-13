import argparse
import os
import shutil
import glob
import ast
import pandas as pd
from utils import *
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



def main():
    parser = argparse.ArgumentParser(description="Check for splice event annotation in CSV.")
    parser.add_argument('-incsv', dest="incsv", help ='Input CSV file path', required=True)
    parser.add_argument('-ref', dest='ref', help="Reference fasta file path", required=True)
    parser.add_argument('-exonlabel', dest='exonlabel', help='Exon relabel CSV file path', required=True)
    parser.add_argument('-fastadir', dest='fastadir', help='out directory  path', required=True)
    args = parser.parse_args()


    if not os.path.exists(args.fastadir):
        os.makedirs(args.fastadir)
        print(f"Directory '{args.fastadir}' was created.")
    else:
        print(f"Directory '{args.fastadir}' already exists.")

    df = pd.read_csv(args.incsv)
    df_exon_unified = pd.read_csv(args.exonlabel)
    fasta_file =  args.ref
    

    try:
        esin = df[df['splice_annotation'] == 'ESin']
        esout = df[df['splice_annotation'] == 'ESout']
        irin = df[df['splice_annotation'] == 'IRin']
        irout = df[df['splice_annotation'] == 'IRout']
        llrs = df[((df['strand'] == '+') & (df['splice_annotation'] == 'A3SSlong')) | 
                     ((df['strand'] == '-') & (df['splice_annotation'] == 'A5SSlong'))]
        lsrsame = df[((df['strand'] == '+') & (df['splice_annotation'] == 'A3SSshort')) | 
                     ((df['strand'] == '-') & (df['splice_annotation'] == 'A5SSshort'))]
        lsrl = df[((df['strand'] == '+') & (df['splice_annotation'] == 'A5SSlong')) | 
                     ((df['strand'] == '-') & (df['splice_annotation'] == 'A3SSlong'))]
        lsamers = df[((df['strand'] == '+') & (df['splice_annotation'] == 'A5SSshort')) | 
                     ((df['strand'] == '-') & (df['splice_annotation'] == 'A3SSshort'))]
    except Exception as e:
        print(f"An error occurred: {e}")

        

    esin_clean = esin_simulation(esin, df_exon_unified)
    esout_clean = esout_simulation(esout, df_exon_unified)
    irin_clean = irin_simulation(irin, df_exon_unified)
    irout_clean = irout_simulation(irout, df_exon_unified)
    llrs_clean = LLRS_simulaion(llrs, df_exon_unified)
    lsrsame_clean = LSRSame_simulation(lsrsame, df_exon_unified)
    lsrl_clean = LSRL_simulation(lsrl, df_exon_unified)
    lsamers_clean = LSameRS_simulation(lsamers, df_exon_unified)
    # true_event_df = pd.concat([esin_clean, esout_clean, irin_clean, irout_clean, llrs_clean, lsamers_clean, lsrl_clean, lsamers_clean])
    # df[df['gene_name'].isin(true_event_df['gene_name'])].to_csv(args.outdir + '/true_event_df.csv', index = False)

    esin_final = esin_true_event(esin_clean, df_exon_unified)
    esout_final = esout_true_event(esout_clean, df_exon_unified)
    irin_final = irin_true_event(irin_clean, df_exon_unified)
    irout_final = irout_true_event(irout_clean, df_exon_unified)
    llrs_final = llrs_true_event(llrs_clean, df_exon_unified)
    lsrsame_final = lsrsame_true_event(lsrsame_clean, df_exon_unified)
    lsrl_final = lsrl_true_event(lsrl_clean, df_exon_unified)
    lsamers_final = lsamers_true_event(lsamers_clean, df_exon_unified)
    true_event_df = pd.concat([esin_final, esout_final, irin_final, irout_final, llrs_final, lsrsame_final, lsrl_final, lsamers_final])
    true_event_df['splice_annotation'] = true_event_df['annotation'] + true_event_df['splice']
    true_event_df.drop_duplicates().to_csv(args.fastadir + '/true_event_df.csv', index = False)

    fasta = pysam.FastaFile(fasta_file)
    esin_ref_fasta, esin_alt_fasta = es_fasta(fasta, esin_clean, df_exon_unified, 'ESin')
    esout_ref_fasta, esout_alt_fasta = es_fasta(fasta, esout_clean, df_exon_unified, 'ESout')
    irin_ref_fasta, irin_alt_fasta = ir_fasta(fasta, irin_clean, df_exon_unified, 'IRin')
    irout_ref_fasta, irout_alt_fasta = ir_fasta(fasta, irout_clean, df_exon_unified, 'IRout')
    LeftLongRightSame_ref_fasta, LeftLongRightSame_alt_fasta = LeftRight_fasta(fasta, llrs_clean, df_exon_unified, 'LeftLongRightSame')
    LeftShortRightSame_ref_fasta, LeftShortRightSame_alt_fasta = LeftRight_fasta(fasta, lsrsame_clean, df_exon_unified, 'LeftShortRightSame')
    LeftSameRightLong_ref_fasta, LeftSameRightLong_alt_fasta = LeftRight_fasta(fasta, lsrl_clean, df_exon_unified, 'LeftSameRightLong')
    LeftSameRightShort_ref_fasta, LeftSameRightShort_alt_fasta = LeftRight_fasta(fasta, lsamers_clean, df_exon_unified, 'LeftSameRightShort')

    ref_concat = esin_ref_fasta + esout_ref_fasta + irin_ref_fasta + irout_ref_fasta + LeftLongRightSame_ref_fasta + LeftShortRightSame_ref_fasta + LeftSameRightLong_ref_fasta + LeftSameRightShort_ref_fasta
    alt_concat = esin_alt_fasta + esout_alt_fasta + irin_alt_fasta + irout_alt_fasta + LeftLongRightSame_alt_fasta + LeftShortRightSame_alt_fasta + LeftSameRightLong_alt_fasta + LeftSameRightShort_alt_fasta

    with open(args.fastadir + '/combined_ref.fasta', "w") as output_handle:
        SeqIO.write(ref_concat, output_handle, "fasta")

    with open(args.fastadir + '/combined_alt.fasta', "w") as output_handle:
        SeqIO.write(alt_concat, output_handle, "fasta")

if __name__ == '__main__':
    main()