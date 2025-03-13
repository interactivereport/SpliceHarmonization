import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import bisect

def esin_simulation(esin, df_exon_unified, N =5):
    esin_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(esin['gene_name'])]
    gene_transcript = esin.merge(esin_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                 right_on = ['Start', 'End', 'gene_name'], 
                                 left_on = ['cassette_2L', 'cassette_2R', 'gene_name'], 
                                 how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]

    splicedin_exon = gene_transcript[['gene_name', 'unified_exon']].drop_duplicates()
    merged_splicedin = esin_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']].merge(splicedin_exon, on=['gene_name', 'unified_exon'])
    merged_splicedin_drop = merged_splicedin.drop_duplicates(['Start', 'End', 'unified_exon'])

    excluded_transcript_id = []
    gene_list = []
    for i, row in merged_splicedin_drop.iterrows():
        splicedin_exon_ = row['unified_exon']
        gene_name = row['gene_name']
        left_ = esin_exon_unified[(esin_exon_unified['gene_name'] == gene_name) & (esin_exon_unified['unified_exon'] <= splicedin_exon_)]
        right_ = esin_exon_unified[(esin_exon_unified['gene_name'] == gene_name) & (esin_exon_unified['unified_exon'] >= splicedin_exon_)]
        if (len(left_) > 0 ) and(len(right_) > 0 ):
            left_filtered = left_[left_['End'] >= left_.iloc[-1]['Start'] ]
            right_filtered = right_[right_['Start'] <=  right_.iloc[0]['End'] ]
            excluded_transcripts_ = list(set(left_filtered['transcript_id']).union(right_filtered['transcript_id']))
            excluded_transcript_id.append(excluded_transcripts_)
            gene_list.append([gene_name]*len(excluded_transcripts_))

    excluded_transcript_id = [item for sublist in excluded_transcript_id for item in sublist]
    gene_list = [item for sublist in gene_list for item in sublist]

    filter_tuples = list(zip(gene_list, excluded_transcript_id))
    filtered_df = esin_exon_unified[~esin_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_tuples)].reset_index(drop=True)

    count_df = filtered_df.groupby(['gene_name', 'transcript_id']).count().rename(columns={'exon_number': 'count'})
    count_df['max_count'] = count_df.groupby('gene_name')['count'].transform(max)
    filtered_df = count_df[count_df['count'] == count_df['max_count']]

    ESin_new = filtered_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('gene_name')
    filtered_ESin = ESin_new.merge(esin_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])

    ESin_sanity = pd.DataFrame(filtered_ESin.groupby('gene_name')['unified_exon'].agg(list)).reset_index().merge(splicedin_exon, left_on ='gene_name', right_on ='gene_name')
    ESin_sanity['unified_label'] = ESin_sanity.apply(lambda x: 'Y' if (x['unified_exon_x'][0] < x['unified_exon_y']) and 
                                                 (x['unified_exon_x'][-1] > x['unified_exon_y']) and (len(x['unified_exon_x']) > N) else 'N', axis=1)
    ESin_final = ESin_sanity[ESin_sanity['unified_label'] == 'Y']

    return ESin_final

def esout_simulation(esout, df_exon_unified, N =5):
    esout_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(esout['gene_name'])]
    gene_transcript = esout.merge(esout_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                  right_on = ['Start', 'End', 'gene_name'], 
                                  left_on = ['cassette_2L', 'cassette_2R', 'gene_name'], 
                                  how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]

    splicedout_exon = gene_transcript[['gene_name', 'unified_exon']].drop_duplicates()


    filter_tuples = list(zip(gene_transcript['gene_name'], gene_transcript['transcript_id']))
    filtered_df = esout_exon_unified[esout_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_tuples)].reset_index(drop=True)

    count_df = filtered_df.groupby(['gene_name', 'transcript_id']).count().rename(columns={'exon_number': 'count'})
    count_df['max_count'] = count_df.groupby('gene_name')['count'].transform(max)
    filtered_df = count_df[count_df['count'] == count_df['max_count']]

    ESout_new = filtered_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('gene_name')
    filtered_ESout = ESout_new.merge(esout_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])

    ESout_sanity = pd.DataFrame(filtered_ESout.groupby('gene_name')['unified_exon'].agg(list)).reset_index().merge(splicedout_exon, left_on ='gene_name', right_on ='gene_name')
    ESout_sanity['unified_label'] = ESout_sanity.apply(lambda x: 'Y' if (x['unified_exon_x'][0] < x['unified_exon_y']) and 
                                                 (x['unified_exon_x'][-1] > x['unified_exon_y']) and (len(x['unified_exon_x']) > N) else 'N', axis=1)
    ESout_final = ESout_sanity[ESout_sanity['unified_label'] == 'Y']

    return ESout_final


def irin_simulation(irin, df_exon_unified, N=5):
    irin_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(irin['gene_name'])]
    gene_transcript_left = irin.merge(irin_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                      right_on = ['Start', 'End', 'gene_name'], 
                                      left_on = ['cassette_1L', 'cassette_1R', 'gene_name'], 
                                      how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript_right = irin.merge(irin_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                       right_on = ['Start', 'End', 'gene_name'], 
                                       left_on = ['cassette_3L', 'cassette_3R', 'gene_name'], 
                                       how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]

    gene_transcript = gene_transcript_left.merge(gene_transcript_right, left_on = ['gene_name', 'transcript_id'], right_on = ['gene_name', 'transcript_id'])
    IRin_exon = gene_transcript[['gene_name', 'unified_exon_x', 'unified_exon_y']].drop_duplicates()
    filter_tuples = list(zip(gene_transcript['gene_name'], gene_transcript['transcript_id']))
    filtered_df = irin_exon_unified[irin_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_tuples)].reset_index(drop=True)

    count_df = filtered_df.groupby(['gene_name', 'transcript_id']).count().rename(columns={'exon_number': 'count'})
    count_df['max_count'] = count_df.groupby('gene_name')['count'].transform(max)
    filtered_df = count_df[count_df['count'] == count_df['max_count']]

    IRin_new = filtered_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('gene_name')
    filtered_IRin = IRin_new.merge(irin_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])

    IRin_sanity = pd.DataFrame(filtered_IRin.groupby('gene_name')['unified_exon'].agg(list)).reset_index().merge(IRin_exon, left_on ='gene_name', right_on ='gene_name')
    IRin_sanity['unified_label'] = IRin_sanity.apply(lambda x: 'Y' if (x['unified_exon'][0] < x['unified_exon_x']) and 
                                                 (x['unified_exon'][-1] > x['unified_exon_y']) and (len(x['unified_exon']) > N) else 'N', axis=1)
    IRin_final = IRin_sanity[IRin_sanity['unified_label'] == 'Y']
    return IRin_final


def irout_simulation(irout, df_exon_unified, N=5):
    irout_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(irout['gene_name'])]
    gene_transcript_left = irout.merge(irout_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                       right_on = ['Start', 'End', 'gene_name'], 
                                       left_on = ['cassette_1L', 'cassette_1R', 'gene_name'], 
                                       how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript_right = irout.merge(irout_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                        right_on = ['Start', 'End', 'gene_name'], 
                                        left_on = ['cassette_3L', 'cassette_3R', 'gene_name'], 
                                        how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript = gene_transcript_left.merge(gene_transcript_right, left_on = ['gene_name', 'transcript_id'], right_on = ['gene_name', 'transcript_id'])
    IRout_exon = gene_transcript[['gene_name', 'unified_exon_x', 'unified_exon_y']].drop_duplicates()

    filter_tuples = list(zip(gene_transcript['gene_name'], gene_transcript['transcript_id']))
    filtered_df = irout_exon_unified[irout_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_tuples)].reset_index(drop=True)
    count_df = filtered_df.groupby(['gene_name', 'transcript_id']).count().rename(columns={'exon_number': 'count'})
    count_df['max_count'] = count_df.groupby('gene_name')['count'].transform(max)
    filtered_df = count_df[count_df['count'] == count_df['max_count']]

    IRout_new = filtered_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('gene_name')
    filtered_IRout = IRout_new.merge(irout_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    IRout_sanity = pd.DataFrame(filtered_IRout.groupby('gene_name')['unified_exon'].agg(list)).reset_index().merge(IRout_exon, left_on ='gene_name', right_on ='gene_name')
    IRout_sanity['unified_label'] = IRout_sanity.apply(lambda x: 'Y' if (x['unified_exon'][0] < x['unified_exon_x']) and 
                                                 (x['unified_exon'][-1] > x['unified_exon_y']) and (len(x['unified_exon']) >= N) else 'N', axis=1)
    
    IRout_final = IRout_sanity[IRout_sanity['unified_label'] == 'Y']
    return IRout_final


def LLRS_simulaion(llrs, df_exon_unified):
    llrs_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(llrs['gene_name'])]
    gene_transcript_ref = llrs.merge(llrs_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                     right_on = ['Start', 'End', 'gene_name'], 
                                     left_on = ['cassette_3L', 'cassette_3R', 'gene_name'], 
                                     how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript_alt = llrs.merge(llrs_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                     right_on = ['Start', 'End', 'gene_name'], 
                                     left_on = ['cassette_2L', 'cassette_2R', 'gene_name'], 
                                     how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    
    unique_ref = set(gene_transcript_ref['transcript_id'].unique())
    unique_alt = set(gene_transcript_alt['transcript_id'].unique())
    exclusive_to_ref = unique_ref - unique_alt
    exclusive_to_alt = unique_alt - unique_ref
    exclusive_ref = gene_transcript_ref[gene_transcript_ref['transcript_id'].isin(exclusive_to_ref)]
    exclusive_alt = gene_transcript_alt[gene_transcript_alt['transcript_id'].isin(exclusive_to_alt)]

    ref_exon = gene_transcript_ref[['gene_name', 'unified_exon']].drop_duplicates()
    alt_exon = gene_transcript_alt[['gene_name', 'unified_exon']].drop_duplicates()

    filter_ref_tuples = list(zip(gene_transcript_ref['gene_name'], gene_transcript_ref['transcript_id']))
    filter_alt_tuples = list(zip(gene_transcript_alt['gene_name'], gene_transcript_alt['transcript_id']))
    filtered_ref_df = llrs_exon_unified[llrs_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_ref_tuples)].reset_index(drop=True)
    filtered_alt_df = llrs_exon_unified[llrs_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_alt_tuples)].reset_index(drop=True)
    llrs_ref_new = filtered_ref_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    llrs_alt_new = filtered_alt_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')

    filtered_ref_llrs = llrs_ref_new.merge(llrs_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    filtered_alt_llrs = llrs_alt_new.merge(llrs_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])

    llrs_sanity_ref = pd.DataFrame(filtered_ref_llrs.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(ref_exon, left_on ='gene_name', right_on ='gene_name')
    llrs_sanity_alt = pd.DataFrame(filtered_alt_llrs.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(alt_exon, left_on ='gene_name', right_on ='gene_name')
    llrs_sanity_ref.columns = ['gene_name', 'transcript_id', 'ref_exon', 'unified_exon_ref']
    llrs_sanity_alt.columns = ['gene_name', 'transcript_id', 'alt_exon', 'unified_exon_alt']

    llrs_sanity_ref['ref_exon_1'] = llrs_sanity_ref.apply(lambda x: sorted(list(set(x['ref_exon']) - set([x['unified_exon_ref']]))), axis=1)
    llrs_sanity_final = llrs_sanity_ref.drop_duplicates('gene_name').merge(llrs_sanity_alt[['gene_name', 'unified_exon_alt']])
    # Ensuring that 'ref_exon_1' contains unique exons and 'unified_exon_alt' is treated as a single-item list before concatenation
    llrs_sanity_final['alt_exon'] = llrs_sanity_final.apply(
        lambda x: sorted(set(x['ref_exon_1']).union([x['unified_exon_alt']])),
        axis=1
    )

    llrs_sanity_final = llrs_sanity_final.drop_duplicates(['gene_name', 'transcript_id'])
    llrs_sanity_final['ref_start'] = llrs_sanity_final.apply(lambda x: llrs_exon_unified[(llrs_exon_unified['gene_name'] == x['gene_name']) & (llrs_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['Start'], axis=1)
    llrs_sanity_final['ref_end'] = llrs_sanity_final.apply(lambda x: llrs_exon_unified[(llrs_exon_unified['gene_name'] == x['gene_name']) & (llrs_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['End'], axis=1)
    llrs_sanity_final['alt_start'] = llrs_sanity_final.apply(lambda x: llrs_exon_unified[(llrs_exon_unified['gene_name'] == x['gene_name']) & (llrs_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['Start'], axis=1)
    llrs_sanity_final['alt_end'] = llrs_sanity_final.apply(lambda x: llrs_exon_unified[(llrs_exon_unified['gene_name'] == x['gene_name']) & (llrs_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['End'], axis=1)

    return llrs_sanity_final


def LSRSame_simulation(lsrsame, df_exon_unified):
    lsrsame_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(lsrsame['gene_name'])]
    gene_transcript_ref = lsrsame.merge(lsrsame_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                        right_on = ['Start', 'End', 'gene_name'], 
                                        left_on = ['cassette_2L', 'cassette_2R', 'gene_name'], 
                                        how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript_alt = lsrsame.merge(lsrsame_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                        right_on = ['Start', 'End', 'gene_name'], 
                                        left_on = ['cassette_3L', 'cassette_3R', 'gene_name'], 
                                        how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    
    unique_ref = set(gene_transcript_ref['transcript_id'].unique())
    unique_alt = set(gene_transcript_alt['transcript_id'].unique())
    exclusive_to_ref = unique_ref - unique_alt
    exclusive_to_alt = unique_alt - unique_ref
    exclusive_ref = gene_transcript_ref[gene_transcript_ref['transcript_id'].isin(exclusive_to_ref)]
    exclusive_alt = gene_transcript_alt[gene_transcript_alt['transcript_id'].isin(exclusive_to_alt)]
    ref_exon = gene_transcript_ref[['gene_name', 'unified_exon']].drop_duplicates()
    alt_exon = gene_transcript_alt[['gene_name', 'unified_exon']].drop_duplicates()
    filter_ref_tuples = list(zip(gene_transcript_ref['gene_name'], gene_transcript_ref['transcript_id']))
    filter_alt_tuples = list(zip(gene_transcript_alt['gene_name'], gene_transcript_alt['transcript_id']))
    filtered_ref_df = lsrsame_exon_unified[lsrsame_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_ref_tuples)].reset_index(drop=True)
    filtered_alt_df = lsrsame_exon_unified[lsrsame_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_alt_tuples)].reset_index(drop=True)
    lsrsame_ref_new = filtered_ref_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    lsrsame_alt_new = filtered_alt_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    filtered_ref_lsrsame = lsrsame_ref_new.merge(lsrsame_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    filtered_alt_lsrsame= lsrsame_alt_new.merge(lsrsame_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    lsrsame_sanity_ref = pd.DataFrame(filtered_ref_lsrsame.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(ref_exon, left_on ='gene_name', right_on ='gene_name')
    lsrsame_sanity_alt = pd.DataFrame(filtered_alt_lsrsame.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(alt_exon, left_on ='gene_name', right_on ='gene_name')
    lsrsame_sanity_ref.columns = ['gene_name', 'transcript_id', 'ref_exon', 'unified_exon_ref']
    lsrsame_sanity_alt.columns = ['gene_name', 'transcript_id', 'alt_exon', 'unified_exon_alt']
    lsrsame_sanity_ref['ref_exon_1'] = lsrsame_sanity_ref.apply(lambda x: sorted(list(set(x['ref_exon']) - set([x['unified_exon_ref']]))), axis=1)
    lsrsame_sanity_final = lsrsame_sanity_ref.drop_duplicates('gene_name').merge(lsrsame_sanity_alt[['gene_name', 'unified_exon_alt']])
    lsrsame_sanity_final['alt_exon'] = lsrsame_sanity_final.apply(
        lambda x: sorted(set(x['ref_exon_1']).union([x['unified_exon_alt']])),
        axis=1
    )

    lsrsame_sanity_final['ref_start'] = lsrsame_sanity_final.apply(lambda x: lsrsame_exon_unified[(lsrsame_exon_unified['gene_name'] == x['gene_name']) & (lsrsame_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['Start'], axis=1)
    lsrsame_sanity_final['ref_end'] = lsrsame_sanity_final.apply(lambda x: lsrsame_exon_unified[(lsrsame_exon_unified['gene_name'] == x['gene_name']) & (lsrsame_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['End'], axis=1)
    lsrsame_sanity_final['alt_start'] = lsrsame_sanity_final.apply(lambda x: lsrsame_exon_unified[(lsrsame_exon_unified['gene_name'] == x['gene_name']) & (lsrsame_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['Start'], axis=1)
    lsrsame_sanity_final['alt_end'] = lsrsame_sanity_final.apply(lambda x: lsrsame_exon_unified[(lsrsame_exon_unified['gene_name'] == x['gene_name']) & (lsrsame_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['End'], axis=1)
    return lsrsame_sanity_final

def LSRL_simulation(lsrl, df_exon_unified):
    lsrl_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(lsrl['gene_name'])]
    gene_transcript_ref = lsrl.merge(lsrl_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                     right_on = ['Start', 'End', 'gene_name'], 
                                     left_on = ['cassette_1L', 'cassette_1R', 'gene_name'], 
                                     how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript_alt = lsrl.merge(lsrl_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], 
                                     right_on = ['Start', 'End', 'gene_name'], 
                                     left_on = ['cassette_2L', 'cassette_2R', 'gene_name'], 
                                     how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    
    unique_ref = set(gene_transcript_ref['transcript_id'].unique())
    unique_alt = set(gene_transcript_alt['transcript_id'].unique())
    exclusive_to_ref = unique_ref - unique_alt
    exclusive_to_alt = unique_alt - unique_ref
    exclusive_ref = gene_transcript_ref[gene_transcript_ref['transcript_id'].isin(exclusive_to_ref)]
    exclusive_alt = gene_transcript_alt[gene_transcript_alt['transcript_id'].isin(exclusive_to_alt)]
    ref_exon = gene_transcript_ref[['gene_name', 'unified_exon']].drop_duplicates()
    alt_exon = gene_transcript_alt[['gene_name', 'unified_exon']].drop_duplicates()
    filter_ref_tuples = list(zip(gene_transcript_ref['gene_name'], gene_transcript_ref['transcript_id']))
    filter_alt_tuples = list(zip(gene_transcript_alt['gene_name'], gene_transcript_alt['transcript_id']))
    filtered_ref_df = lsrl_exon_unified[lsrl_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_ref_tuples)].reset_index(drop=True)
    filtered_alt_df = lsrl_exon_unified[lsrl_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_alt_tuples)].reset_index(drop=True)
    lsrl_ref_new = filtered_ref_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    lsrl_alt_new = filtered_alt_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    filtered_ref_lsrl = lsrl_ref_new.merge(lsrl_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    filtered_alt_lsrl = lsrl_alt_new.merge(lsrl_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    lsrl_sanity_ref = pd.DataFrame(filtered_ref_lsrl.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(ref_exon, left_on ='gene_name', right_on ='gene_name')
    lsrl_sanity_alt = pd.DataFrame(filtered_alt_lsrl.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(alt_exon, left_on ='gene_name', right_on ='gene_name')
    lsrl_sanity_ref.columns = ['gene_name', 'transcript_id', 'ref_exon', 'unified_exon_ref']
    lsrl_sanity_alt.columns = ['gene_name', 'transcript_id', 'alt_exon', 'unified_exon_alt']
    lsrl_sanity_ref['ref_exon_1'] = lsrl_sanity_ref.apply(lambda x: sorted(list(set(x['ref_exon']) - set([x['unified_exon_ref']]))), axis=1)
    lsrl_sanity_final = lsrl_sanity_ref.drop_duplicates('gene_name').merge(lsrl_sanity_alt[['gene_name', 'unified_exon_alt']])
    # Ensuring that 'ref_exon_1' contains unique exons and 'unified_exon_alt' is treated as a single-item list before concatenation
    lsrl_sanity_final['alt_exon'] = lsrl_sanity_final.apply(
        lambda x: sorted(set(x['ref_exon_1']).union([x['unified_exon_alt']])),
        axis=1
    )

    lsrl_sanity_final = lsrl_sanity_final.drop_duplicates(['gene_name', 'transcript_id'])
    lsrl_sanity_final['ref_start'] = lsrl_sanity_final.apply(lambda x: lsrl_exon_unified[(lsrl_exon_unified['gene_name'] == x['gene_name']) & (lsrl_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['Start'], axis=1)
    lsrl_sanity_final['ref_end'] = lsrl_sanity_final.apply(lambda x: lsrl_exon_unified[(lsrl_exon_unified['gene_name'] == x['gene_name']) & (lsrl_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['End'], axis=1)
    lsrl_sanity_final['alt_start'] = lsrl_sanity_final.apply(lambda x: lsrl_exon_unified[(lsrl_exon_unified['gene_name'] == x['gene_name']) & (lsrl_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['Start'], axis=1)
    lsrl_sanity_final['alt_end'] = lsrl_sanity_final.apply(lambda x: lsrl_exon_unified[(lsrl_exon_unified['gene_name'] == x['gene_name']) & (lsrl_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['End'], axis=1)

    return lsrl_sanity_final

def LSameRS_simulation(lsamers, df_exon_unified):
    lsamers_exon_unified = df_exon_unified[df_exon_unified['gene_name'].isin(lsamers['gene_name'])]
    gene_transcript_ref = lsamers.merge(lsamers_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], right_on = ['Start', 'End', 'gene_name'], left_on = ['cassette_2L', 'cassette_2R', 'gene_name'], how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    gene_transcript_alt = lsamers.merge(lsamers_exon_unified[['Start', 'End', 'gene_name', 'transcript_id', 'unified_exon']], right_on = ['Start', 'End', 'gene_name'], left_on = ['cassette_1L', 'cassette_1R', 'gene_name'], how = 'inner')[['gene_name', 'transcript_id', 'unified_exon']]
    unique_ref = set(gene_transcript_ref['transcript_id'].unique())
    unique_alt = set(gene_transcript_alt['transcript_id'].unique())
    exclusive_to_ref = unique_ref - unique_alt
    exclusive_to_alt = unique_alt - unique_ref
    exclusive_ref = gene_transcript_ref[gene_transcript_ref['transcript_id'].isin(exclusive_to_ref)]
    exclusive_alt = gene_transcript_alt[gene_transcript_alt['transcript_id'].isin(exclusive_to_alt)]
    ref_exon = gene_transcript_ref[['gene_name', 'unified_exon']].drop_duplicates()
    alt_exon = gene_transcript_alt[['gene_name', 'unified_exon']].drop_duplicates()
    filter_ref_tuples = list(zip(gene_transcript_ref['gene_name'], gene_transcript_ref['transcript_id']))
    filter_alt_tuples = list(zip(gene_transcript_alt['gene_name'], gene_transcript_alt['transcript_id']))
    filtered_ref_df = lsamers_exon_unified[lsamers_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_ref_tuples)].reset_index(drop=True)
    filtered_alt_df = lsamers_exon_unified[lsamers_exon_unified.set_index(['gene_name', 'transcript_id']).index.isin(filter_alt_tuples)].reset_index(drop=True)
    lsamers_ref_new = filtered_ref_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    lsamers_alt_new = filtered_alt_df.reset_index()[['gene_name', 'transcript_id']].drop_duplicates('transcript_id')
    filtered_ref_lsmaers = lsamers_ref_new.merge(lsamers_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    filtered_alt_lsmaers = lsamers_alt_new.merge(lsamers_exon_unified[['gene_name', 'transcript_id', 'unified_exon']])
    lsamers_sanity_ref = pd.DataFrame(filtered_ref_lsmaers.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(ref_exon, left_on ='gene_name', right_on ='gene_name')
    lsamers_sanity_alt = pd.DataFrame(filtered_alt_lsmaers.groupby(['gene_name', 'transcript_id'])['unified_exon'].agg(list)).reset_index().merge(alt_exon, left_on ='gene_name', right_on ='gene_name')
    lsamers_sanity_ref.columns = ['gene_name', 'transcript_id', 'ref_exon', 'unified_exon_ref']
    lsamers_sanity_alt.columns = ['gene_name', 'transcript_id', 'alt_exon', 'unified_exon_alt']
    lsamers_sanity_ref['ref_exon_1'] = lsamers_sanity_ref.apply(lambda x: sorted(list(set(x['ref_exon']) - set([x['unified_exon_ref']]))), axis=1)

    lsamers_sanity_final = lsamers_sanity_ref.drop_duplicates('gene_name').merge(lsamers_sanity_alt[['gene_name', 'unified_exon_alt']])
    # Ensuring that 'ref_exon_1' contains unique exons and 'unified_exon_alt' is treated as a single-item list before concatenation
    lsamers_sanity_final['alt_exon'] = lsamers_sanity_final.apply(
        lambda x: sorted(set(x['ref_exon_1']).union([x['unified_exon_alt']])),
        axis=1
    )

    lsamers_sanity_final = lsamers_sanity_final.drop_duplicates(['gene_name', 'transcript_id'])
    lsamers_sanity_final['ref_start'] = lsamers_sanity_final.apply(lambda x: lsamers_exon_unified[(lsamers_exon_unified['gene_name'] == x['gene_name']) & (lsamers_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['Start'], axis=1)
    lsamers_sanity_final['ref_end'] = lsamers_sanity_final.apply(lambda x: lsamers_exon_unified[(lsamers_exon_unified['gene_name'] == x['gene_name']) & (lsamers_exon_unified['unified_exon'] == x['unified_exon_ref']) ].iloc[0]['End'], axis=1)
    lsamers_sanity_final['alt_start'] = lsamers_sanity_final.apply(lambda x: lsamers_exon_unified[(lsamers_exon_unified['gene_name'] == x['gene_name']) & (lsamers_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['Start'], axis=1)
    lsamers_sanity_final['alt_end'] = lsamers_sanity_final.apply(lambda x: lsamers_exon_unified[(lsamers_exon_unified['gene_name'] == x['gene_name']) & (lsamers_exon_unified['unified_exon'] == x['unified_exon_alt']) ].iloc[0]['End'], axis=1)

    return lsamers_sanity_final



def es_fasta(fasta, es, df_exon_unified, annotation):
    records_ref = []
    for i, row in es.iterrows():
        gene_name = row['gene_name']
        exons = row['unified_exon_x']
        sub_ = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'].isin(exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
        chr_ = sub_.iloc[0]['Chromosome']
        strand = sub_.iloc[0]['Strand']
        seq_ = ''
    
        tst = list(zip(sub_['Start'] -1, sub_['End']))
        if strand == '-':
            tst = tst[::-1]
        for (start, end) in tst:
            sequence = fasta.fetch(chr_, start, end)
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            seq_ += sequence
        
        seq_record = SeqRecord(
                    Seq(seq_),
                    id=f"{gene_name}:{sub_['Start'].min()}:{sub_['End'].max()}:{strand}:{annotation}",
                )
        records_ref.append(seq_record)


    records_alt = []
    if annotation == 'ESin':
        for i, row in es.iterrows():
            gene_name = row['gene_name']
            exons = sorted(row['unified_exon_x'] + [row['unified_exon_y']])
            sub_ = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                        (df_exon_unified['unified_exon'].isin(exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
            chr_ = sub_.iloc[0]['Chromosome']
            strand = sub_.iloc[0]['Strand']
            seq_ = ''
            tst = list(zip(sub_['Start'] -1, sub_['End']))
            if strand == '-':
                tst = tst[::-1]
            for (start, end) in tst:
                sequence = fasta.fetch(chr_, start, end)
                if strand == '-':
                    sequence = str(Seq(sequence).reverse_complement())
                seq_ += sequence
            
            seq_record = SeqRecord(
                        Seq(seq_),
                        id=f"{gene_name}:{sub_['Start'].min()}:{sub_['End'].max()}:{strand}:{annotation}",
                    )
            records_alt.append(seq_record)


    elif annotation == 'ESout':
        for i, row in es.iterrows():
            gene_name = row['gene_name']
            exons = sorted(list(set(row['unified_exon_x']) - set([row['unified_exon_y']])))
            sub_ = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                        (df_exon_unified['unified_exon'].isin(exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
            chr_ = sub_.iloc[0]['Chromosome']
            strand = sub_.iloc[0]['Strand']
            seq_ = ''
            tst = list(zip(sub_['Start'] -1, sub_['End']))
            if strand == '-':
                tst = tst[::-1]
            for (start, end) in tst:
                sequence = fasta.fetch(chr_, start, end)
                if strand == '-':
                    sequence = str(Seq(sequence).reverse_complement())
                seq_ += sequence
            
            seq_record = SeqRecord(
                        Seq(seq_),
                        id=f"{gene_name}:{sub_['Start'].min()}:{sub_['End'].max()}:{strand}:{annotation}",
                    )
            records_alt.append(seq_record)

    else:
        print('No ES splice events')

    return records_ref, records_alt
    


def ir_fasta(fasta, ir, df_exon_unified, annotation):
    records = []
    for i, row in ir.iterrows():
        gene_name = row['gene_name']
        exons = row['unified_exon']
        sub_ = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'].isin(exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
        chr_ = sub_.iloc[0]['Chromosome']
        strand = sub_.iloc[0]['Strand']
        seq_ = ''
    
        tst = list(zip(sub_['Start'] -1, sub_['End']))
        if strand == '-':
            tst = tst[::-1]
        for (start, end) in tst:
            sequence = fasta.fetch(chr_, start, end)
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            seq_ += sequence
        
        seq_record = SeqRecord(
                    Seq(seq_),
                    id=f"{gene_name}:{sub_['Start'].min()}:{sub_['End'].max()}:{strand}:{annotation}",
                )
        records.append(seq_record)

    if annotation == 'IRin':
        irin_records_ref = records
    elif annotation == 'IRout':
        irout_records_alt = records
    else:
        print('Not IR splice events')

    records= []
    for i, row in ir.iterrows():
        gene_name = row['gene_name']
        left_exons = row['unified_exon'][:row['unified_exon'].index(row['unified_exon_x'])]
        right_exons = row['unified_exon'][row['unified_exon'].index(row['unified_exon_y']) + 1 :]
        
        sub_left = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'].isin(left_exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
        sub_right = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'].isin(right_exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
        
        sub_ = pd.DataFrame({'gene_name': [gene_name], 'Start': [df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'] == row['unified_exon_x'])]['Start'].iloc[0]], 'End' : [df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'] == row['unified_exon_y'])]['End'].iloc[0]]})
        chr_ = sub_left.iloc[0]['Chromosome']
        strand = sub_left.iloc[0]['Strand']
        seq_ = ''
        tst_left = list(zip(sub_left['Start'] -1, sub_left['End']))
        tst_ =  list(zip(sub_['Start'] -1, sub_['End']))
        tst_right= list(zip(sub_right['Start'] -1, sub_right['End']))
        concatenated_list = tst_left + tst_ + tst_right
        
        if strand == '-':
            concatenated_list = concatenated_list[::-1]
        for (start, end) in concatenated_list:
            sequence = fasta.fetch(chr_, start, end)
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            seq_ += sequence
        
        seq_record = SeqRecord(
                    Seq(seq_),
                    id=f"{gene_name}:{sub_left['Start'].min()}:{sub_right['End'].max()}:{strand}:{annotation}",
                )
        records.append(seq_record)


    if annotation == 'IRin':
        irin_records_alt = records
        return irin_records_ref, irin_records_alt
    elif annotation == 'IRout':
        irout_records_ref = records
        return irout_records_ref, irout_records_alt
    else:
        print('Not IR splice events')


def LeftRight_fasta(fasta, df, df_exon_unified, annotation):
    records_ref = []
    for i, row in df.iterrows():
        gene_name = row['gene_name']
        exons = row['ref_exon']
        sub_ = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'].isin(exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
        chr_ = sub_.iloc[0]['Chromosome']
        strand = sub_.iloc[0]['Strand']
        seq_ = ''
    
        tst = list(zip(sub_['Start'] -1, sub_['End']))
        if strand == '-':
            tst = tst[::-1]
        for (start, end) in tst:
            sequence = fasta.fetch(chr_, start, end)
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            seq_ += sequence
        
        seq_record = SeqRecord(
                    Seq(seq_),
                    id=f"{gene_name}:{sub_['Start'].min()}:{sub_['End'].max()}:{strand}:{annotation}",
                )
        records_ref.append(seq_record)

    records_alt = []
    for i, row in df.iterrows():
        gene_name = row['gene_name']
        exons = row['alt_exon']
        sub_ = df_exon_unified[(df_exon_unified['gene_name'] == gene_name) &
                    (df_exon_unified['unified_exon'].isin(exons))][['Chromosome', 'Strand', 'Start', 'End', 'unified_exon']].drop_duplicates()
        chr_ = sub_.iloc[0]['Chromosome']
        strand = sub_.iloc[0]['Strand']
        seq_ = ''
        tst = list(zip(sub_['Start'] -1, sub_['End']))
        if strand == '-':
            tst = tst[::-1]
        for (start, end) in tst:
            sequence = fasta.fetch(chr_, start, end)
            if strand == '-':
                sequence = str(Seq(sequence).reverse_complement())
            seq_ += sequence
        
        seq_record = SeqRecord(
                    Seq(seq_),
                    id=f"{gene_name}:{sub_['Start'].min()}:{sub_['End'].max()}:{strand}:{annotation}",
                )
        records_alt.append(seq_record)

    return records_ref, records_alt



    

def leafcutter_base_metrics(intron_file, true_event_file, cutoffs):
    true_event_df = pd.read_csv(true_event_file)
    df = pd.read_csv(intron_file)
    df = df[df['gene'] != '.']
    
    # first pass - with R is the same, L is + 1
    leafcutter_res = []
    for idx, group in df.groupby('clusterID'):
        clustername, genename = group['clusterID'].unique()[0], group['gene'].unique()[0]
        current_event = true_event_df[true_event_df['gene_name'] == group['gene'].iloc[0]]
        if len(current_event) > 0 :
            if current_event['annotation'].iloc[0] == 'ES':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0],current_event['cassette_2R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                            'End': [current_event['cassette_2L'].iloc[0] + 1 ,current_event['cassette_3L'].iloc[0] + 1 , current_event['cassette_3L'].iloc[0] + 1] })
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                if len(group_merged) >= 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'ES', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0]  + 1 , current_event['cassette_3L'].iloc[0] + 1]})

                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len A5SS + is {}'.format(len(group_merged)))
                print(genename)
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0] + 1 , current_event['cassette_3L'].iloc[0] + 1]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A5SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0]  + 1 , current_event['cassette_3L'].iloc[0] + 1]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0] + 1 , current_event['cassette_3L'].iloc[0] + 1 ]})  
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS + is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif  current_event['annotation'].iloc[0] == 'IR':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0] + 1 ]}) 
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len IR + is {}'.format(len(group_merged)))
                if len(group_merged) == 1:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'IR', leafcutter_dpsi])   
                
            else:
                # leafcutter_res.append([clustername, genename, 'NAN', float('nan')])  
                leafcutter_res.append([clustername, genename, 'NAN', group_merged['deltapsi'].abs().max()]) 
                
    leafcutter_with_true_p1 = pd.DataFrame(leafcutter_res, columns=['clusterID', 'gene', 'true_annotation', 'dpsi_max'])

    # first pass - with R is the same, L is the same as the +1bp buffering room
    leafcutter_res = []
    for idx, group in df.groupby('clusterID'):
        clustername, genename = group['clusterID'].unique()[0], group['gene'].unique()[0]
        current_event = true_event_df[true_event_df['gene_name'] == group['gene'].iloc[0]]
        if len(current_event) > 0 :
            if current_event['annotation'].iloc[0] == 'ES':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0],current_event['cassette_2R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                            'End': [current_event['cassette_2L'].iloc[0],current_event['cassette_3L'].iloc[0] , current_event['cassette_3L'].iloc[0] ] })
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                if len(group_merged) >= 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'ES', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0], current_event['cassette_3L'].iloc[0]]})

                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len A5SS + is {}'.format(len(group_merged)))
                # print(genename)
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0]  , current_event['cassette_3L'].iloc[0] ]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A5SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0]  , current_event['cassette_3L'].iloc[0] ]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0] , current_event['cassette_3L'].iloc[0]]})  
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS + is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif  current_event['annotation'].iloc[0] == 'IR':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0] + 1 ]}) 
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len IR + is {}'.format(len(group_merged)))
                if len(group_merged) == 1:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'IR', leafcutter_dpsi])   
                
            else:
                # leafcutter_res.append([clustername, genename, 'NAN', float('nan')]) 
                leafcutter_res.append([clustername, genename, 'NAN', group_merged['deltapsi'].abs().max()])  
                
    leafcutter_with_true_p2 = pd.DataFrame(leafcutter_res, columns=['clusterID', 'gene', 'true_annotation', 'dpsi_max'])

    leafcutter_with_true_all = pd.concat([ leafcutter_with_true_p1, leafcutter_with_true_p2])
    leafcutter_with_true_all = leafcutter_with_true_all.sort_values(by='dpsi_max', ascending=False).drop_duplicates(subset=[ 'gene', 'true_annotation'], keep='first')
    leafcutter_with_true_all['dpsi_max'] = leafcutter_with_true_all['dpsi_max'].fillna(-1)
    leafcutter_with_true_merged = true_event_df.merge(leafcutter_with_true_all, left_on = ['gene_name'], right_on = ['gene'], how = 'outer')
    # leafcutter_with_true_merged['strand_new'] = leafcutter_with_true_merged['clusterID'].apply(lambda x: x.split('_')[-1] if isinstance(x, str) else np.nan)

    leafcutter_base_tpr, leafcutter_base_fpr, leafcutter_base_precision, leafcutter_base_recall, leafcutter_base_f1 = [], [], [], [], []
    leafcutter_base_counts = []
    for cutoff in cutoffs:
        leafcutter_base_tp_counts, leafcutter_base_tn_counts, leafcutter_base_fp_counts, leafcutter_base_fn_counts = \
        len(leafcutter_with_true_merged[(~leafcutter_with_true_merged['gene'].isna()) & 
        (leafcutter_with_true_merged['dpsi_max'] >= cutoff) 
            & (~leafcutter_with_true_merged['gene_name'].isna()) ]), \
        len(leafcutter_with_true_merged[(leafcutter_with_true_merged['gene_name'].isna()) & ((leafcutter_with_true_merged['gene'].isna()) | 
                                                (leafcutter_with_true_merged['dpsi_max'] < cutoff))]), \
        len(leafcutter_with_true_merged[((leafcutter_with_true_merged['gene_name'].isna()) | (leafcutter_with_true_merged['true_annotation'] == 'NAN')) & 
                                                                    (leafcutter_with_true_merged['dpsi_max']>= cutoff )]) , \
        len(leafcutter_with_true_merged[((leafcutter_with_true_merged['gene'].isna()) | 
                                                (leafcutter_with_true_merged['dpsi_max'] < cutoff)) & (~leafcutter_with_true_merged['gene_name'].isna()) ])
        # print(leafcutter_base_tp_counts+leafcutter_base_tn_counts+leafcutter_base_fp_counts+leafcutter_base_fn_counts)
        leafcutter_base_fp_counts = max(0, df[df['deltapsi'] >=  cutoff]['clusterID'].nunique() - leafcutter_base_tp_counts)
        leafcutter_base_tpr.append((leafcutter_base_tp_counts)/(leafcutter_base_tp_counts+leafcutter_base_fn_counts))
        leafcutter_base_fpr.append((leafcutter_base_fp_counts)/(leafcutter_base_fp_counts+leafcutter_base_tn_counts) if (leafcutter_base_fp_counts+leafcutter_base_tn_counts)!= 0 else float('nan'))
        leafcutter_base_precision.append((leafcutter_base_tp_counts)/(leafcutter_base_tp_counts+leafcutter_base_fp_counts) if (leafcutter_base_tp_counts+leafcutter_base_fp_counts)!=0 else float('nan'))
        leafcutter_base_recall.append((leafcutter_base_tp_counts)/(leafcutter_base_tp_counts+leafcutter_base_fn_counts))
        leafcutter_base_f1.append(2*leafcutter_base_precision[-1]*leafcutter_base_recall[-1]/(leafcutter_base_precision[-1] + leafcutter_base_recall[-1]) if (leafcutter_base_precision[-1] + leafcutter_base_recall[-1])!=0 else float('nan'))
        leafcutter_base_counts.append([leafcutter_base_tp_counts, leafcutter_base_tn_counts, leafcutter_base_fp_counts, leafcutter_base_fn_counts])


    return leafcutter_base_tpr, leafcutter_base_fpr, leafcutter_base_precision, leafcutter_base_recall, leafcutter_base_f1, leafcutter_base_counts

def leafcutter_base_truedf(intron_file, true_event_file, cutoff):
    true_event_df = pd.read_csv(true_event_file)
    df = pd.read_csv(intron_file)
    df = df[df['gene'] != '.']
    
    # first pass - with R is the same, L is + 1
    leafcutter_res = []
    for idx, group in df.groupby('clusterID'):
        clustername, genename = group['clusterID'].unique()[0], group['gene'].unique()[0]
        current_event = true_event_df[true_event_df['gene_name'] == group['gene'].iloc[0]]
        if len(current_event) > 0 :
            if current_event['annotation'].iloc[0] == 'ES':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0],current_event['cassette_2R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                            'End': [current_event['cassette_2L'].iloc[0] + 1 ,current_event['cassette_3L'].iloc[0] + 1 , current_event['cassette_3L'].iloc[0] + 1] })
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                if len(group_merged) >= 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'ES', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0]  + 1 , current_event['cassette_3L'].iloc[0] + 1]})

                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len A5SS + is {}'.format(len(group_merged)))
                # print(genename)
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0] + 1 , current_event['cassette_3L'].iloc[0] + 1]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A5SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0]  + 1 , current_event['cassette_3L'].iloc[0] + 1]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0] + 1 , current_event['cassette_3L'].iloc[0] + 1 ]})  
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS + is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif  current_event['annotation'].iloc[0] == 'IR':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0] + 1 ]}) 
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len IR + is {}'.format(len(group_merged)))
                if len(group_merged) == 1:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'IR', leafcutter_dpsi])   
                
            else:
                # leafcutter_res.append([clustername, genename, 'NAN', float('nan')])  
                leafcutter_res.append([clustername, genename, 'NAN', group_merged['deltapsi'].abs().max()])
                
    leafcutter_with_true_p1 = pd.DataFrame(leafcutter_res, columns=['clusterID', 'gene', 'true_annotation', 'dpsi_max'])

    # first pass - with R is the same, L is the same as the +1bp buffering room
    leafcutter_res = []
    for idx, group in df.groupby('clusterID'):
        clustername, genename = group['clusterID'].unique()[0], group['gene'].unique()[0]
        current_event = true_event_df[true_event_df['gene_name'] == group['gene'].iloc[0]]
        if len(current_event) > 0 :
            if current_event['annotation'].iloc[0] == 'ES':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0],current_event['cassette_2R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                            'End': [current_event['cassette_2L'].iloc[0],current_event['cassette_3L'].iloc[0] , current_event['cassette_3L'].iloc[0] ] })
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                if len(group_merged) >= 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'ES', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0], current_event['cassette_3L'].iloc[0]]})

                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len A5SS + is {}'.format(len(group_merged)))
                # print(genename)
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A5SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0]  , current_event['cassette_3L'].iloc[0] ]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A5SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A5SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '-'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_2R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0]  , current_event['cassette_3L'].iloc[0] ]})
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS - is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif (current_event['annotation'].iloc[0] == 'A3SS') and (current_event['strand'].iloc[0] == '+'):
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0], current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_2L'].iloc[0] , current_event['cassette_3L'].iloc[0]]})  
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])

                # print('len A3SS + is {}'.format(len(group_merged)))
                if len(group_merged) == 2:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'A3SS', leafcutter_dpsi])

            elif  current_event['annotation'].iloc[0] == 'IR':
                current_event_intron = pd.DataFrame({'Start': [current_event['cassette_1R'].iloc[0]], 
                                                    'End': [current_event['cassette_3L'].iloc[0] + 1 ]}) 
                group_merged = group.merge(current_event_intron, left_on = ['start', 'end'], right_on = ['Start', 'End'])
                # print('len IR + is {}'.format(len(group_merged)))
                if len(group_merged) == 1:
                    leafcutter_dpsi = group_merged['deltapsi'].abs().max()
                else:
                    leafcutter_dpsi = float('nan')
                leafcutter_res.append([clustername, genename, 'IR', leafcutter_dpsi])   
                
            else:
                # leafcutter_res.append([clustername, genename, 'NAN', float('nan')])  
                leafcutter_res.append([clustername, genename, 'NAN', group_merged['deltapsi'].abs().max()])
                
    leafcutter_with_true_p2 = pd.DataFrame(leafcutter_res, columns=['clusterID', 'gene', 'true_annotation', 'dpsi_max'])

    leafcutter_with_true_all = pd.concat([ leafcutter_with_true_p1, leafcutter_with_true_p2])
    leafcutter_with_true_all = leafcutter_with_true_all.sort_values(by='dpsi_max', ascending=False).drop_duplicates(subset=[ 'gene', 'true_annotation'], keep='first')
    leafcutter_with_true_all['dpsi_max'] = leafcutter_with_true_all['dpsi_max'].fillna(-1)
    leafcutter_with_true_merged = true_event_df.merge(leafcutter_with_true_all, left_on = ['gene_name'], right_on = ['gene'], how = 'outer')
    # leafcutter_with_true_merged['strand_new'] = leafcutter_with_true_merged['clusterID'].apply(lambda x: x.split('_')[-1] if isinstance(x, str) else np.nan)
    leafcutter_base_tp= \
    leafcutter_with_true_merged[(~leafcutter_with_true_merged['gene'].isna()) & 
        (leafcutter_with_true_merged['dpsi_max'] >= cutoff) 
            & (~leafcutter_with_true_merged['gene_name'].isna()) ]


    return leafcutter_base_tp
    
def find_bounds(sorted_list, x):
    """
    Given a sorted list and a value x, return the last element in the list
    that is strictly less than x (left) and the first element that is strictly greater than x (right).
    If no such element exists, return None for that bound.
    """
    # # Find the insertion index for x in the sorted list.
    left_index = bisect.bisect_left(sorted_list, x) - 1  # Last index where element < x
    right_index = bisect.bisect_right(sorted_list, x)      # First index where element > x

    left_value = sorted_list[left_index] if left_index >= 0 else None
    right_value = sorted_list[right_index] if right_index < len(sorted_list) else None

    return pd.Series({'left': left_value, 'right': right_value})


def esin_true_event(esin_clean, df_exon_unified):
    esin_clean[['left', 'right']] = esin_clean.apply(lambda row: find_bounds(row['unified_exon_x'], row['unified_exon_y']), axis=1)
    esin_cassette2 = esin_clean[['gene_name', 'unified_exon_y']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'unified_exon_y'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    esin_cassette2 = esin_cassette2.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_2L', 'End':'cassette_2R'})
    esin_cassette1 = esin_clean[['gene_name', 'left']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'left'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    esin_cassette1 = esin_cassette1.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_1L', 'End':'cassette_1R'})
    esin_cassette3 = esin_clean[['gene_name', 'right']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'right'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    esin_cassette3 = esin_cassette3.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_3L', 'End':'cassette_3R'})
    esin_final = esin_cassette1.merge(esin_cassette2, on = ['gene_name', 'chr', 'strand']).merge(esin_cassette3,on = ['gene_name', 'chr', 'strand'])[['gene_name', 'chr', 'strand',
                                                                                                                                        'cassette_1L', 'cassette_1R', 
                                                                                                                                        'cassette_2L', 'cassette_2R',
                                                                                                                                        'cassette_3L', 'cassette_3R', ]]
    esin_final['annotation'] = 'ES'
    esin_final['splice'] = 'in'

    return esin_final.drop_duplicates()

def esout_true_event(esout_clean, df_exon_unified):
    esout_clean[['left', 'right']] = esout_clean.apply(lambda row: find_bounds(row['unified_exon_x'], row['unified_exon_y']), axis=1)
    esout_cassette2 = esout_clean[['gene_name', 'unified_exon_y']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'unified_exon_y'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    esout_cassette2 = esout_cassette2.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_2L', 'End':'cassette_2R'})
    esout_cassette1 = esout_clean[['gene_name', 'left']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'left'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    esout_cassette1 = esout_cassette1.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_1L', 'End':'cassette_1R'})
    esout_cassette3 = esout_clean[['gene_name', 'right']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'right'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    esout_cassette3 = esout_cassette3.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_3L', 'End':'cassette_3R'})

    esout_final = esout_cassette1.merge(esout_cassette2, on = ['gene_name', 'chr', 'strand']).merge(esout_cassette3,on = ['gene_name', 'chr', 'strand'])[['gene_name', 'chr', 'strand',
                                                                                                                                        'cassette_1L', 'cassette_1R', 
                                                                                                                                        'cassette_2L', 'cassette_2R',
                                                                                                                                        'cassette_3L', 'cassette_3R', ]]
    esout_final['annotation'] = 'ES'
    esout_final['splice'] = 'out'
    return esout_final.drop_duplicates()


def irin_true_event(irin_clean, df_exon_unified):
    irin_cassette1 = irin_clean[['gene_name', 'unified_exon_x']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                  left_on = ['gene_name', 'unified_exon_x'],
                                                 right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    irin_cassette1 = irin_cassette1.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_1L', 'End':'cassette_1R'})

    irin_cassette3 = irin_clean[['gene_name', 'unified_exon_y']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'unified_exon_y'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    irin_cassette3 = irin_cassette3.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_3L', 'End':'cassette_3R'})

    irin_final = irin_cassette1.merge(irin_cassette3, on = ['gene_name', 'chr', 'strand'])[['gene_name', 'chr', 'strand','cassette_1L', 'cassette_1R', 'cassette_3L', 'cassette_3R', ]]
    irin_final['cassette_2L'] = irin_final['cassette_1L'] 
    irin_final['cassette_2R'] = irin_final['cassette_3R'] 
    irin_final['annotation'] = 'IR'
    irin_final['splice'] = 'in'
    return irin_final.drop_duplicates()

def irout_true_event(irout_clean, df_exon_unified):
    irout_cassette1 = irout_clean[['gene_name', 'unified_exon_x']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                  left_on = ['gene_name', 'unified_exon_x'],
                                                 right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    irout_cassette1 = irout_cassette1.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_1L', 'End':'cassette_1R'})

    irout_cassette3 = irout_clean[['gene_name', 'unified_exon_y']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'unified_exon_y'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    irout_cassette3 = irout_cassette3.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_3L', 'End':'cassette_3R'})

    irout_final = irout_cassette1.merge(irout_cassette3, on = ['gene_name', 'chr', 'strand'])[['gene_name', 'chr', 'strand','cassette_1L', 'cassette_1R', 'cassette_3L', 'cassette_3R', ]]
    irout_final['cassette_2L'] = irout_final['cassette_1L'] 
    irout_final['cassette_2R'] = irout_final['cassette_3R'] 
    irout_final['annotation'] = 'IR'
    irout_final['splice'] = 'out'
    return irout_final.drop_duplicates()


def llrs_true_event(llrs_clean, df_exon_unified):
    llrs_clean[['left', 'right']] = llrs_clean.apply(lambda row: find_bounds(row['ref_exon'], row['unified_exon_ref']), axis=1)
    llrs_cassette1 = llrs_clean[['gene_name', 'left']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'left'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    llrs_cassette1 = llrs_cassette1.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_1L', 'End':'cassette_1R'})
    llrs_clean['cassette_2L'] = llrs_clean['alt_start']
    llrs_clean['cassette_3L'] = llrs_clean['ref_start']
    llrs_clean['cassette_2R'] =  llrs_clean['cassette_3R']  = llrs_clean['ref_end']
    llrs_final = llrs_cassette1.merge(llrs_clean, on = ['gene_name'])[['gene_name', 'chr', 'strand','cassette_1L', 'cassette_1R', 'cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R', ]]
    llrs_final['splice'] = 'long'
    llrs_final['annotation'] = llrs_final['strand'].apply(lambda x: 'A3SS' if x == '+' else 'A5SS')

    return llrs_final.drop_duplicates()

def lsrsame_true_event(lsrsame_clean, df_exon_unified):
    lsrsame_clean[['left', 'right']] = lsrsame_clean.apply(lambda row: find_bounds(row['ref_exon'], row['unified_exon_ref']), axis=1)
    lsrsame_cassette1 = lsrsame_clean[['gene_name', 'left']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'left'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    lsrsame_cassette1 = lsrsame_cassette1.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_1L', 'End':'cassette_1R'})
    lsrsame_clean['cassette_2L'] = lsrsame_clean['ref_start']
    lsrsame_clean['cassette_3L'] = lsrsame_clean['alt_start']
    lsrsame_clean['cassette_2R'] =  lsrsame_clean['cassette_3R']  = lsrsame_clean['ref_end']
    lsrsame_final = lsrsame_cassette1.merge(lsrsame_clean, on = ['gene_name'])[['gene_name', 'chr', 'strand','cassette_1L', 'cassette_1R', 'cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R', ]]
    lsrsame_final['splice'] = 'short'
    lsrsame_final['annotation'] = lsrsame_final['strand'].apply(lambda x: 'A3SS' if x == '+' else 'A5SS')
    return lsrsame_final.drop_duplicates()


def lsrl_true_event(lsrl_clean, df_exon_unified):
    lsrl_clean[['left', 'right']] = lsrl_clean.apply(lambda row: find_bounds(row['ref_exon'], row['unified_exon_ref']), axis=1)
    lsrl_cassette3= lsrl_clean[['gene_name', 'right']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'right'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    lsrl_cassette3 = lsrl_cassette3.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_3L', 'End':'cassette_3R'})
    lsrl_clean['cassette_1R'] = lsrl_clean['ref_end']
    lsrl_clean['cassette_2R'] = lsrl_clean['alt_end']
    lsrl_clean['cassette_1L'] =  lsrl_clean['cassette_2L']  = lsrl_clean['ref_start']
    lsrl_final = lsrl_cassette3.merge(lsrl_clean, on = ['gene_name'])[['gene_name', 'chr', 'strand','cassette_1L', 'cassette_1R', 'cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R', ]]
    lsrl_final['splice'] = 'long'
    lsrl_final['annotation'] = lsrl_final['strand'].apply(lambda x: 'A5SS' if x == '+' else 'A3SS')
    return lsrl_final.drop_duplicates()


def lsamers_true_event(lsamers_clean, df_exon_unified):
    lsamers_clean[['left', 'right']] = lsamers_clean.apply(lambda row: find_bounds(row['ref_exon'], row['unified_exon_ref']), axis=1)
    lsamers_cassette3= lsamers_clean[['gene_name', 'right']].merge(df_exon_unified[['Start', 'End', 'gene_name', 'unified_exon', 'Chromosome', 'Strand']].drop_duplicates(), 
                                                    left_on = ['gene_name', 'right'],
                                                    right_on = ['gene_name', 'unified_exon'], how ='left')[['gene_name',  'Start', 'End', 'unified_exon','Chromosome', 'Strand']]
    lsamers_cassette3 = lsamers_cassette3.rename(columns={'Chromosome':'chr', 'Strand':'strand', 'Start': 'cassette_3L', 'End':'cassette_3R'})
    lsamers_clean['cassette_1R'] = lsamers_clean['alt_end']
    lsamers_clean['cassette_2R'] = lsamers_clean['ref_end']
    lsamers_clean['cassette_1L'] =  lsamers_clean['cassette_2L']  = lsamers_clean['ref_start']
    lsamers_final = lsamers_cassette3.merge(lsamers_clean, on = ['gene_name'])[['gene_name', 'chr', 'strand','cassette_1L', 'cassette_1R', 'cassette_2L', 'cassette_2R', 'cassette_3L', 'cassette_3R', ]]
    lsamers_final['splice'] = 'short'
    lsamers_final['annotation'] = lsamers_final['strand'].apply(lambda x: 'A5SS' if x == '+' else 'A3SS')
    return lsamers_final.drop_duplicates()





def majiq_true_ES(ES, true_event_df):
    ES['cassette_list'] = ES['lsv_id'].apply(lambda x: [sub.split(':')[-1] for sub in x.split(";")])
    ES['noswap'] = ES['cassette_list'].apply(lambda x: int(x[0].split("-")[0]) < int(x[1].split("-")[0]))
    ES['1L'] = ES[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][0].split("-")[0]) if x['noswap'] else int(x['cassette_list'][1].split("-")[0]), axis=1)
    ES['1R'] = ES[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][0].split("-")[1]) if x['noswap'] else int(x['cassette_list'][1].split("-")[1]), axis=1)
    ES['3L'] = ES[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][1].split("-")[0]) if x['noswap'] else int(x['cassette_list'][0].split("-")[0]), axis=1)
    ES['3R'] = ES[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][1].split("-")[1]) if x['noswap'] else int(x['cassette_list'][0].split("-")[1]), axis=1)
    ES_with_prediction_all = true_event_df[true_event_df['annotation'] == 'ES'].merge(ES, on = ['gene_name', 'strand'], how = 'outer')
    es_fp1 = ES_with_prediction_all[ES_with_prediction_all['cassette_1L'].isna()]
    es_fn1 = ES_with_prediction_all[ES_with_prediction_all['1L'].isna()]
    ES_with_prediction_inter = true_event_df[true_event_df['annotation'] == 'ES'].merge(ES, on = ['gene_name', 'strand'], how = 'inner')
    ES_with_prediction_inter['d_diff_sum']  = ES_with_prediction_inter.apply(lambda x: np.abs(x['cassette_1L'] - x['1L']) + np.abs(x['cassette_1R'] - x['1R']) +  np.abs(x['cassette_3L'] - x['3L']) + np.abs(x['cassette_3R'] - x['3R']), axis=1)
    es_tp = ES_with_prediction_inter[ES_with_prediction_inter['d_diff_sum'] <= 12]
    es_fp2 = ES_with_prediction_inter[ES_with_prediction_inter['d_diff_sum'] > 12]
    es_fn2 = ES_with_prediction_inter[ES_with_prediction_inter['d_diff_sum'] > 12]
    tp_counts, tn_counts, fn_counts, fp_counts = len(es_tp), 0, len(es_fn1) + len(es_fn2), len(es_fp1) + len(es_fp2)
    return tp_counts,  tn_counts, fn_counts, fp_counts, es_tp

def majiq_true_ASS(ASS, true_event_df):
    ASS['cassette_list'] = ASS['lsv_id'].apply(lambda x: [sub.split(':')[-1] for sub in x.split(";")])
    ASS['noswap'] = ASS['cassette_list'].apply(lambda x: int(x[0].split("-")[0]) < int(x[1].split("-")[0]) if len(x) > 1 else 'NAN')
    ASS['1L'] = ASS[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][0].split("-")[0]) if x['noswap'] else int(x['cassette_list'][1].split("-")[0]), axis=1)
    ASS['1R'] = ASS[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][0].split("-")[1]) if x['noswap'] else int(x['cassette_list'][1].split("-")[1]), axis=1)
    ASS['3L'] = ASS[['cassette_list', 'noswap']].apply(
        lambda x: int(x['cassette_list'][1].split("-")[0]) 
            if x['noswap'] == True 
            else int(x['cassette_list'][0].split("-")[0]) 
            if x['noswap'] != 'NAN' 
            else 'NAN',
        axis=1
    )
    ASS['3R'] = ASS[['cassette_list', 'noswap']].apply(
        lambda x: int(x['cassette_list'][1].split("-")[1]) 
            if x['noswap'] == True 
            else int(x['cassette_list'][0].split("-")[1]) 
            if x['noswap'] != 'NAN' 
            else 'NAN',
        axis=1
    )
    ass_with_prediction_all = true_event_df[true_event_df['annotation'].isin(['A5SS', 'A3SS'])].merge(ASS, on = ['gene_name', 'strand'], how = 'outer')
    ass_fp1 = ass_with_prediction_all[ass_with_prediction_all['cassette_1L'].isna()]
    ass_fn1 = ass_with_prediction_all[ass_with_prediction_all['1L'].isna()]
    ass_with_prediction_inter = true_event_df[true_event_df['annotation'].isin(['A5SS', 'A3SS'])].merge(ASS, on = ['gene_name', 'strand'], how = 'inner')
    ass_with_prediction_inter['d_diff_sum']  = ass_with_prediction_inter.apply(lambda x: min([np.abs(x['cassette_2L'] - x['1L']) + np.abs(x['cassette_2R'] - x['1R']), (np.abs(x['cassette_1L'] - x['1L']) + np.abs(x['cassette_1R'] - x['1R'])), (np.abs(x['cassette_3L'] - x['1L']) + np.abs(x['cassette_3R'] - x['1R'])),
                                                                                                (np.abs(x['cassette_2L'] - x['3L']) + np.abs(x['cassette_2R'] - x['3R'])),(np.abs(x['cassette_1L'] - x['3L']) + np.abs(x['cassette_1R'] - x['3R'])), (np.abs(x['cassette_3L'] - x['3L']) + np.abs(x['cassette_3R'] - x['3R']))]) if x['noswap'] != 'NAN'
                                                                                                else min([(np.abs(x['cassette_2L'] - x['1L']) + np.abs(x['cassette_2R'] - x['1R'])), (np.abs(x['cassette_1L'] - x['1L']) + np.abs(x['cassette_1R'] - x['1R'])), (np.abs(x['cassette_3L'] - x['1L']) + np.abs(x['cassette_3R'] - x['1R']))]) , axis=1)
    ass_tp = ass_with_prediction_inter[ass_with_prediction_inter['d_diff_sum'] <= 8]
    ass_fp2 = ass_with_prediction_inter[ass_with_prediction_inter['d_diff_sum'] > 8]
    ass_fn2 = ass_with_prediction_inter[ass_with_prediction_inter['d_diff_sum'] > 8]
    tp_counts, tn_counts, fn_counts, fp_counts = len(ass_tp), 0, len(ass_fn1) + len(ass_fn2), len(ass_fp1) + len(ass_fp2)
    return tp_counts, tn_counts, fn_counts, fp_counts, ass_tp

def majiq_true_IR(IR, true_event_df):
    IR['cassette_list'] = IR['lsv_id'].apply(lambda x: [sub.split(':')[-1] for sub in x.split(";")])
    IR['noswap'] = IR['cassette_list'].apply(lambda x: int(x[0].split("-")[0]) < int(x[1].split("-")[0]) if len(x) > 1 else 'NAN')
    IR['1L'] = IR[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][0].split("-")[0]) if x['noswap'] else int(x['cassette_list'][1].split("-")[0]), axis=1)
    IR['1R'] = IR[['cassette_list', 'noswap']].apply(lambda x: int(x['cassette_list'][0].split("-")[1]) if x['noswap'] else int(x['cassette_list'][1].split("-")[1]), axis=1)
    ir_with_prediction_all = true_event_df[true_event_df['annotation'] == 'IR'].merge(IR, on = ['gene_name', 'strand'], how = 'outer')
    ir_fp1 = ir_with_prediction_all[ir_with_prediction_all['cassette_1L'].isna()]
    ir_fn1 = ir_with_prediction_all[ir_with_prediction_all['1L'].isna()]
    ir_with_prediction_inter = true_event_df[true_event_df['annotation'] == 'IR'].merge(IR, on = ['gene_name', 'strand'], how = 'inner')
    ir_with_prediction_inter['d_diff_sum']  = ir_with_prediction_inter.apply(lambda x: min([np.abs(x['cassette_2L'] - x['1L']) + np.abs(x['cassette_2R'] - x['1R']), (np.abs(x['cassette_1L'] - x['1L']) + np.abs(x['cassette_1R'] - x['1R'])), (np.abs(x['cassette_3L'] - x['1L']) + np.abs(x['cassette_3R'] - x['1R'])),
                                                                                                (np.abs(x['cassette_2L'] - x['3L']) + np.abs(x['cassette_2R'] - x['3R'])),(np.abs(x['cassette_1L'] - x['3L']) + np.abs(x['cassette_1R'] - x['3R'])), (np.abs(x['cassette_3L'] - x['3L']) + np.abs(x['cassette_3R'] - x['3R']))]) if x['noswap'] != 'NAN'
                                                                                                else min([(np.abs(x['cassette_2L'] - x['1L']) + np.abs(x['cassette_2R'] - x['1R'])), (np.abs(x['cassette_1L'] - x['1L']) + np.abs(x['cassette_1R'] - x['1R'])), (np.abs(x['cassette_3L'] - x['1L']) + np.abs(x['cassette_3R'] - x['1R']))]) , axis=1)
    ir_tp = ir_with_prediction_inter[ir_with_prediction_inter['d_diff_sum'] <= 4]
    ir_fp2 = ir_with_prediction_inter[ir_with_prediction_inter['d_diff_sum'] > 4]
    ir_fn2 = ir_with_prediction_inter[ir_with_prediction_inter['d_diff_sum'] > 4]
    tp_counts, tn_counts, fn_counts, fp_counts = len(ir_tp), 0, len(ir_fn1) + len(ir_fn2), len(ir_fp1) + len(ir_fp2)
    return tp_counts, tn_counts, fn_counts, fp_counts, ir_tp



def majiq_baseline_true_event(true_event_df, df):
    ES = df[~df['cassette'].isna()]
    A3SS = df[~df['alt3'].isna()]
    A5SS = df[~df['alt5'].isna()]
    IR = df[~df['ir'].isna()]

    es_tp, es_tn, es_fn, es_fp  = 0,0,0,0
    a5ss_tp, a5ss_tn, a5ss_fn, a5ss_fp =0,0,0, 0
    a3ss_tp, a3ss_tn, a3ss_fn, a3ss_fp =0, 0, 0,0
    ir_tp, ir_tn, ir_fn, ir_fp = 0, 0, 0,0

    if len(ES) > 0:
        es_tp, es_tn, es_fn, es_fp, es_tp_events = majiq_true_ES(ES, true_event_df)
    if len(A5SS) > 0:
        a5ss_tp, a5ss_tn, a5ss_fn, a5ss_fp , a5ss_tp_events= majiq_true_ASS(A5SS, true_event_df)
    if len(A3SS) > 0:
        a3ss_tp, a3ss_tn, a3ss_fn, a3ss_fp, a3ss_tp_events = majiq_true_ASS(A3SS, true_event_df)
    if len(IR) > 0: 
        ir_tp, ir_tn, ir_fn, ir_fp, ir_tp_events = majiq_true_IR(IR, true_event_df)

    tp_counts = es_tp +a5ss_tp + a3ss_tp + ir_tp
    tn_counts  = 0
    fn_counts = es_fn + a5ss_fn + a3ss_fn + ir_fn
    fp_counts = df['num_events'].sum() - tp_counts

    return tp_counts, tn_counts, fn_counts, fp_counts


def majiq_baseline_TP_res(true_event_df, df):
    ES = df[~df['cassette'].isna()]
    A3SS = df[~df['alt3'].isna()]
    A5SS = df[~df['alt5'].isna()]
    IR = df[~df['ir'].isna()]

    es_tp, es_tn, es_fn, es_fp  = 0,0,0,0
    a5ss_tp, a5ss_tn, a5ss_fn, a5ss_fp =0,0,0, 0
    a3ss_tp, a3ss_tn, a3ss_fn, a3ss_fp =0, 0, 0,0
    ir_tp, ir_tn, ir_fn, ir_fp = 0, 0, 0,0

    if len(ES) > 0:
        es_tp, es_tn, es_fn, es_fp, es_tp_events = majiq_true_ES(ES, true_event_df)
    if len(A5SS) > 0:
        a5ss_tp, a5ss_tn, a5ss_fn, a5ss_fp , a5ss_tp_events= majiq_true_ASS(A5SS, true_event_df)
    if len(A3SS) > 0:
        a3ss_tp, a3ss_tn, a3ss_fn, a3ss_fp, a3ss_tp_events = majiq_true_ASS(A3SS, true_event_df)
    if len(IR) > 0: 
        ir_tp, ir_tn, ir_fn, ir_fp, ir_tp_events = majiq_true_IR(IR, true_event_df)

    return es_tp_events, a5ss_tp_events, a3ss_tp_events, ir_tp_events


