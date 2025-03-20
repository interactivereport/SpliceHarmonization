import numpy as np
import itertools
import pandas as pd
import warnings
import re
from itertools import combinations
import ast
import glob 
import pysam
from collections import Counter
from scipy.stats import ttest_ind, combine_pvalues
#from statsmodels.stats.multitest import multipletests
from concurrent.futures import ProcessPoolExecutor
import json
from multiprocessing import Pool, cpu_count

def left_right_combine(left_list, right_list):
    # es_list[es] = [coverage_list, left_list, right_list]
    listR_dict = {key: [value for _, value in group] for key, group in itertools.groupby(right_list, key=lambda x: x[0])}
    
    output = [subleft_list + [value] for subleft_list in left_list for value in listR_dict.get(subleft_list[-1], [])]
    
    output = set(tuple(x) for x in output)
    if output:
        return list(output)
    else:
        return None
    
def create_graphs(connections):
    # Initialize dictionary to store neighbors
    graph_dict = {}
    
    # Populate the dictionary with connections
    for connection in connections:
        node1, node2 = connection
        graph_dict.setdefault(node1, []).append(node2)
        graph_dict.setdefault(node2, []).append(node1)
    
    # Initialize set for visited nodes
    visited = set()
    
    # List to store the resulting graphs
    graphs = []
    
    # Perform depth-first search
    def dfs(node, graph):
        graph.append(node)
        visited.add(node)
        neighbors = graph_dict.get(node, [])
        for neighbor in neighbors:
            if neighbor not in visited:
                dfs(neighbor, graph)
    
    # Iterate through each node in the dictionary
    for node in graph_dict:
        if node not in visited:
            # Create an empty graph
            graph = []
            # Perform DFS starting from the current node
            dfs(node, graph)
            # Append the graph to the list of graphs
            graphs.append(graph)
    
    return graphs

def find_triplets_graph(linked_graph):
    if len(linked_graph) < 3:
        return [[int(x) for x in linked_graph]]
    triplets = [list(int(x) for x in triplet) for triplet in combinations(linked_graph, 3)]
    triples_index = list(combinations(range(len(linked_graph)), 3))
    #return triplets, triples_index
    return triplets

def nodes_tothree(node_list):
    node_dict = {}
    for entry in node_list:
        key, value = entry[0], entry[1]
        if key not in node_dict:
            node_dict[key] = set()
        node_dict[key].add(value)

    return node_dict

def create_threenode_graph(connections):
    def dfs(connections, start, depth, visited=None, path=None):
        if visited is None:
            visited = set()
        if path is None:
            path = [start]
        else:
            path.append(start)
        visited.add(start)
        if depth == 0:
            return [path]
        paths = []
        for connection in connections:
            node1, node2 = connection
            if node1 == start and node2 not in visited:
                sub_paths = dfs(connections, node2, depth - 1, visited.copy(), path[:])
                paths.extend(sub_paths)
            elif node2 == start and node1 not in visited:
                sub_paths = dfs(connections, node1, depth - 1, visited.copy(), path[:])
                paths.extend(sub_paths)
        return paths
    
    start_node_candidate  = list(set([x for sublist in connections for x in sublist]))
    result  = [] 
    for start_node in start_node_candidate:
        res_ = dfs(connections , start_node, 2)
        for re in res_:
            if re and sorted(re) not in result:
                result.append(sorted(re))
    return result

def dfs_backup(graph, start, depth, visited=None):
    if visited is None:
        visited = set()
    visited.add(start)
    if depth == 0:
        return [[start]]
    paths = []
    for neighbor in graph.get(start, []):
        if neighbor not in visited:
            sub_paths = dfs_backup(graph, neighbor, depth - 1, visited.copy())
            for path in sub_paths:
                paths.append([start] + path)
    visited.remove(start)
    return paths

def exonlist_to_threenode_graph(exon_list, main_graph):
    exon_list =  [[int(element) for element in list(y)] for y in set([tuple(x) for x in exon_list])]
    graph = [int(x) for x in main_graph]
    graph_dict = nodes_tothree(exon_list)
    #graph_dict = create_threenode_graph(exon_list)
    
    nodes_to_delete = []
    for node in graph_dict:
        if node not in graph:
            nodes_to_delete.append(node)

    for node in nodes_to_delete:
        del graph_dict[node]

    start_node = sorted(graph_dict.keys())[0]
    three_node_graph = list(itertools.chain.from_iterable([dfs_backup(graph_dict, start_node, 2) for start_node in sorted(graph_dict.keys())]))
    
    if not three_node_graph:
        pairwise_pairs = [(key, value) for key, values in graph_dict.items() for value in values]
        pairwise_three_node_graph = create_threenode_graph(pairwise_pairs)
        if not pairwise_three_node_graph:
            return pairwise_pairs
        else:
            return pairwise_three_node_graph
        # #return create_threenode_graph(pairwise_pairs)
        # return pairwise_pairs
    else:
        return three_node_graph

class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, x, y):
        root_x = self.find(x)
        root_y = self.find(y)
        if root_x != root_y:
            self.parent[root_x] = root_y


def get_counts(bam_file_paths, threads=8):
                with Pool(threads) as pool:
                                results = pool.map(bam_count, bam_file_paths)
                                # assumes each file path is unique in the list
                                results = {result[1]: result[0] for result in results}
                return results  


def bam_count(bam_file_path):
                with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
                                # return a tuple for traceability in async processes
                                return (round(bam_file.count()/1e6), bam_file_path)


### without Parallel
# def process_bam_count(bam_file_paths):
#     alignment_count = []
#     for bam_file_path in bam_file_paths:
#         with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
#             alignment_count.append(round(bam_file.count()/1e6) ) 
#     return alignment_countstatsmodels

def count_alignments(bam_file_path):
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        return round(bam_file.count() / 1e6)

def process_bam_count(bam_file_paths):
    with Pool() as pool:
        alignment_counts = pool.map(count_alignments, bam_file_paths)
    return alignment_counts


def process_read(read, position, row, junction_dict):
    """Helper function to process each read and update junction counts."""
    for op, length in read.cigar:
        if op in {0, 2}:  # M (match/mismatch) or D (deletion) operations
            position += length
        elif op == 3:  # N (skipped region from splicing)
            junction_start = position
            junction_end = position + length
            # Check specific junctions based on the row and update counts
            check_and_update_junction(junction_start, junction_end, row, junction_dict)
            position += length

def check_and_update_junction(junction_start, junction_end, row, junction_dict):
    """Check if junction matches any predefined cassettes and update counts."""
    cassettes = [(row['cassette_1R'], row['cassette_3L'], '13'),
                 (row['cassette_1R'], row['cassette_2L'], '12'),
                 (row['cassette_2R'], row['cassette_3L'], '23')]
    for start, end, key in cassettes:
        if junction_start == start and junction_end == end:
            junction_key = (junction_start, junction_end)
            if junction_key in junction_dict[key]:
                junction_dict[key][junction_key] += 1
            else:
                junction_dict[key][junction_key] = 1

def median_of_three(a, b, c):
    return sorted([a, b, c])[1]


# def junctioncounts_df_new(df, alt_bam_file_paths, ref_bam_file_paths):
#     junction_counts = {'13': {}, '12': {}, '23': {}}
#     for bam_file_path in alt_bam_file_paths:
#         with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
#             for idx, row in df.iterrows():
#                 for read in bam_file.fetch(row['chr'], row['cassette_1R'], row['cassette_3L']):
#                     if read.is_unmapped or not any(op == 3 for op, length in read.cigar):
#                         continue
#                     process_read(read, read.reference_start, row, junction_counts)

#     df['alt_count13'] = df.apply(lambda x: junction_counts['13'].get((x['cassette_1R'], x['cassette_3L']), 0), axis=1)
#     df['alt_count12'] = df.apply(lambda x: junction_counts['12'].get((x['cassette_1R'], x['cassette_2L']), 0), axis=1)
#     df['alt_count23'] = df.apply(lambda x: junction_counts['23'].get((x['cassette_2R'], x['cassette_3L']), 0), axis=1)

#     for bam_file_path in ref_bam_file_paths:
#         with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
#             for idx, row in df.iterrows():
#                 for read in bam_file.fetch(row['chr'], row['cassette_1R'], row['cassette_3L']):
#                     if read.is_unmapped or not any(op == 3 for op, length in read.cigar):
#                         continue
#                     process_read(read, read.reference_start, row, junction_counts)

#     df['ref_count13'] = df.apply(lambda x: junction_counts['13'].get((x['cassette_1R'], x['cassette_3L']), 0), axis=1)
#     df['ref_count12'] = df.apply(lambda x: junction_counts['12'].get((x['cassette_1R'], x['cassette_2L']), 0), axis=1)
#     df['ref_count23'] = df.apply(lambda x: junction_counts['23'].get((x['cassette_2R'], x['cassette_3L']), 0), axis=1)


#     return df

def process_bam_files(df, bam_file_paths, column_suffix, depth_counts):
    junction_counts = {'13': {}, '12': {}, '23': {}}

    for file_idx, bam_file_path in enumerate(bam_file_paths):
        alignment_count = depth_counts[file_idx]
        with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
            for idx, row in df.iterrows():
                try:
                    for read in bam_file.fetch(row['chr'], row['cassette_1R'], row['cassette_3L']):
                        if read.is_unmapped or not any(op == 3 for op, length in read.cigar):
                                continue
                        process_read(read, read.reference_start, row, junction_counts)
                except ValueError:
                    # Handle possible ValueError if fetch parameters are incorrect
                    continue

            # for key in junction_counts:
            #     for position, count in junction_counts[key].items():
            #         ratio = count / alignment_count if alignment_count > 0 else 0
            #         if position in junction_ratios[key]:
            #             junction_ratios[key][position] += ratio
            #         else:
            #             junction_ratios[key][position] = ratio

    # Update DataFrame directly with counts from junction_counts
    for key in junction_counts:
        df[f'{column_suffix}_count{key}'] = df.apply(lambda x: junction_counts[key].get((x[f'cassette_{key[0]}R'], x[f'cassette_{key[1]}L']), 0)/np.sum(depth_counts)/len(depth_counts), axis=1)

### without parallel
# def process_bam_files_countTolist(df, bam_file_paths, column_suffix):
#     for key in ['13', '12', '23']:
#         df[f'{column_suffix}_count{key}'] = [[] for _ in range(len(df))]
#     for bam_file_path in bam_file_paths:
#         junction_counts = {'13': {}, '12': {}, '23': {}}
#         with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
#             for idx, row in df.iterrows():
#                 try:
#                     for read in bam_file.fetch(row['chr'], row['cassette_1R'], row['cassette_3L']):
#                         if read.is_unmapped or not any(op == 3 for op, length in read.cigar):
#                                 continue
#                         process_read(read, read.reference_start, row, junction_counts)
#                 except ValueError:
#                     # Handle possible ValueError if fetch parameters are incorrect
#                     continue    
#         for key in junction_counts:
#             results = df.apply(lambda x: junction_counts[key].get((x[f'cassette_{key[0]}R'], x[f'cassette_{key[1]}L']), 0), axis=1)
#             for i, result in enumerate(results):
#                 df.at[i, f'{column_suffix}_count{key}'].append(result)


def process_single_bam(bam_file_path, df_subset):
    junction_counts = {'13': {}, '12': {}, '23': {}}
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for index, row in df_subset.iterrows():
            try:
                for read in bam_file.fetch(row['chr'], row['cassette_1R'], row['cassette_3L']):
                    if read.is_unmapped or not any(op == 3 for op, length in read.cigar):
                        continue
                    process_read(read, read.reference_start, row, junction_counts)
            except ValueError:
                continue
    return junction_counts

def process_bam_files_countTolist(df, bam_file_paths, column_suffix):
    # Initialize empty lists in the dataframe for storing counts
    for key in ['13', '12', '23']:
        df[f'{column_suffix}_count{key}'] = [[] for _ in range(len(df))]

    # Prepare subsets of the dataframe to pass to the pool
    subsets = [(bam_file_path, df.copy()) for bam_file_path in bam_file_paths]

    with Pool() as pool:
        results = pool.starmap(process_single_bam, subsets)

    # Aggregate results
    for junction_counts in results:
        for key in ['13', '12', '23']:
            for i, row in df.iterrows():
                count = junction_counts[key].get((row[f'cassette_{key[0]}R'], row[f'cassette_{key[1]}L']), 0)
                df.at[i, f'{column_suffix}_count{key}'].append(count)

def junctioncounts_df_new(df, alt_bam_file_paths, ref_bam_file_paths, alt_depth, ref_depth):

    process_bam_files(df, alt_bam_file_paths, 'alt', alt_depth)
    process_bam_files(df, ref_bam_file_paths, 'ref', ref_depth)
    return df

def junctioncounts_df_countTolist(df, alt_bam_file_paths, ref_bam_file_paths):
    process_bam_files_countTolist(df, alt_bam_file_paths, 'alt')
    process_bam_files_countTolist(df, ref_bam_file_paths, 'ref')
    return df
# def junctioncounts_df(df, bam_file_paths):

#     ref_junction_counts13 = {}
#     ref_junction_counts12 = {}
#     ref_junction_counts23 = {}

#     for bam_file_path in bam_file_paths:
#         bam_file = pysam.AlignmentFile(bam_file_path, "rb")
#         for idx, row in df.iterrows():
#             for read in bam_file.fetch(row['chr'], row['cassette_1R'], row['cassette_3L']):
#                 if read.is_unmapped or not any((op == 3) for op, length in read.cigar):
#                     continue  # Skip unmapped reads and reads without a 'N' operation
#                 position = read.reference_start  # Start position of the read
#                 for op, length in read.cigar:
#                     if op == 0 or op == 2:  # M or D operations, move along the reference
#                         position += length
#                     elif op == 3:  # N operation, indicates a junction
#                         junction_start = position
#                         junction_end = position + length

#                         if (junction_start ==  row['cassette_1R'] ) or (junction_end == row['cassette_3L'] ):
#                             junction_key = (junction_start, junction_end)
#                         if junction_key in ref_junction_counts13:
#                             ref_junction_counts13[junction_key] += 1
#                         else:
#                             ref_junction_counts13[junction_key] = 1
                            
#                         if (junction_start ==  row['cassette_1R'] ) or (junction_end == row['cassette_2L'] ):
#                             junction_key = (junction_start, junction_end)
#                             if junction_key in ref_junction_counts12:
#                                 ref_junction_counts12[junction_key] += 1
#                             else:
#                                 ref_junction_counts12[junction_key] = 1
                                
#                         if (junction_start ==  row['cassette_2R'] ) or (junction_end == row['cassette_3L'] ):
#                             junction_key = (junction_start, junction_end)
#                             if junction_key in ref_junction_counts23:
#                                 ref_junction_counts23[junction_key] += 1
#                             else:
#                                 ref_junction_counts23[junction_key] = 1

#                         position += length
#             bam_file.close()

#     df['count13'] = df.apply(lambda x: ref_junction_counts13.get((x['cassette_1R'], x['cassette_3L']), 0), axis=1)
#     df['count12'] = df.apply(lambda x: ref_junction_counts12.get((x['cassette_1R'], x['cassette_2L']), 0), axis=1)
#     df['count23'] = df.apply(lambda x: ref_junction_counts23.get((x['cassette_2R'], x['cassette_3L']), 0), axis=1)

#     return df


def junctionfilter(df, median_cutoff):

    ir_df = df[df['type'] == 'IR']

    # Process 'ES' type rows, calculate the median, and filter
    es_df = df[df['type'] == 'ES']
    es_df['counts_sum'] = es_df.apply(lambda row: median_of_three(row['count13'], row['count12'], row['count23']), axis=1)
    es_df = es_df[es_df['counts_sum'] > median_cutoff]
    es_df = es_df.drop(['counts_sum'], axis=1)

    # Process rows that are not 'IR' or 'ES', find unique counts and their minimum
    other_df = df[~df['type'].isin(['IR', 'ES'])]
    other_df['min_unique_count'] = other_df.apply(lambda row: min(set([row['count1'], row['count2'], row['count3']])), axis=1)
    other_df = other_df[other_df['min_unique_count'] > median_cutoff]
    other_df = other_df.drop(['min_unique_count'], axis=1)

    # Combine all the filtered DataFrames
    final_df = pd.concat([ir_df, es_df, other_df])

    return final_df



# Function to load and parse the JSON file
def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)  # Load and parse the JSON content into a dictionary
    return data


def perform_t_test(control, treatment):
    stat, p_value = ttest_ind(control, treatment, equal_var= False)
    return p_value

# def perform_t_test(data):
#     # Unpack data
#     ref_count13, alt_count13, ref_count12, alt_count12, ref_count23, alt_count23 = data

#     # Ensure data is in the correct format (list) before performing t-test
#     def safe_eval(data_item):
#         if isinstance(data_item, str):
#             return eval(data_item)
#         else:
#             return data_item
    
#     ref_count13 = safe_eval(ref_count13)
#     alt_count13 = safe_eval(alt_count13)
#     ref_count12 = safe_eval(ref_count12)
#     alt_count12 = safe_eval(alt_count12)
#     ref_count23 = safe_eval(ref_count23)
#     alt_count23 = safe_eval(alt_count23)

#     # Perform t-tests
#     pvalue13 = ttest_ind(np.array(ref_count13), np.array(alt_count13), equal_var=False).pvalue
#     pvalue12 = ttest_ind(np.array(ref_count12), np.array(alt_count12), equal_var=False).pvalue
#     pvalue23 = ttest_ind(np.array(ref_count23), np.array(alt_count23), equal_var=False).pvalue
#     return [pvalue13, pvalue12, pvalue23]

# def perform_t_test_strtolist(control, treatment):
#     stat, p_value = ttest_ind(ast.literal_eval(control), ast.literal_eval(treatment), equal_var= False)
#     return p_value

# def convert_str_to_list(s):
#     return ast.literal_eval(s)


def FisherP(multiple_p_values, method='fisher'):
    combined_stat, combined_p_value = combine_pvalues(multiple_p_values, method=method)
    return combined_p_value 

def pairwise_division(row, normalize_counts):
    return [x / y for x, y in zip(row, normalize_counts)]

def pairwise_division_new(row, normalize_counts):
    # Ensure both row and normalize_counts contain numeric values
    row = [float(x) for x in row]
    normalize_counts = [float(y) for y in normalize_counts]
    return [x / y for x, y in zip(row, normalize_counts)]

def remove_nans(lst):
    return [x for x in lst if pd.notna(x)]


def max_in_row(row):
    combined_list = sum((row[col] for col in row.index), [])
    return max(combined_list)


def F1_label_prep(df, dpsi_cutoff, score_cutoff, fisher_pvalue, fc_cutoff, maxcount_cutoff, inputjson):

    fc_thred, invfc_thred = fc_cutoff, 1/fc_cutoff
    libdepths = load_json(inputjson)

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

    # data_junctioncount = [(row['ref_count13'], row['alt_count13'], row['ref_count12'], row['alt_count12'], row['ref_count23'], row['alt_count23']) 
    #     for index, row in df.iterrows()]


    # with ProcessPoolExecutor() as executor:
    #     res_= list(executor.map(perform_t_test, data_junctioncount))
    
    # df['ttest_pvalue'] = res_
    # df['ttest_pvalue_removeNAs'] = df['ttest_pvalue'].apply(remove_nans)
    # df['fisher_pvalue'] = df['ttest_pvalue_removeNAs'].apply(FisherP)
    #df['count_normalize'] = df[['alt_count13', 'alt_count12', 'alt_count23', 'ref_count13', 'ref_count12', 'ref_count23']].max(axis=1)
    # df['count_normalize'] = df[['alt_count13', 'alt_count12', 'alt_count23', 'ref_count13', 'ref_count12', 'ref_count23']].apply(max_in_row, axis=1)

    val = maxcount_cutoff/np.mean(list(libdepths['alt_depth'].values()) + list(libdepths['ref_depth'].values()))
    df['count_normalize_label'] = df['count_normalize'] >= val


    df['junctionDE_label'] = False
    df.loc[(df['fisher_pvalue'] <= fisher_pvalue) & 
                          (df['count_normalize_label'] == True) & 
                          (df['junctionLabel'] == True),'junctionDE_label'] = True
    df['psi_label'] = False
    df.loc[(df['score'] > score_cutoff) & (df['dpsi'].abs() > dpsi_cutoff), 'psi_label'] = True
    df['F1_label'] = df.apply( lambda x: 'TP' if (x['junctionDE_label'] and x['psi_label'] ) else
                                                            'TN' if (not x['junctionDE_label'] and not x['psi_label']) else
                                                            'FP' if (not x['junctionDE_label'] and x['psi_label']) else
                                                            'FN' if (x['junctionDE_label'] and not x['psi_label']) else 
                                                            'NA',
                                                            axis=1)
    rmats = df[df['method'] == 'rmats']
    leafcutter = df[df['method'] == 'leafcutter']
    majiq = df[df['method'] == 'majiq']

    return pd.concat([rmats,leafcutter, majiq ])


def unique_event_output(df):
    result = df.groupby(['event_label']).agg({
        'method':  lambda x: list(set(x)),
        'score': 'mean',
        'dpsi':'mean',
        'fisher_pvalue': 'mean',
        'junctionLabel': lambda x: Counter(list(x)).most_common(1)[0][0],
        'count_normalize_label': lambda x: list(set(x))[0],
    }).reset_index()
    
    return result

def F1_label(df, dpsi_cutoff, score_cutoff, fisher_pvalue):

    result = df.groupby(['event_label']).agg({
    'method':  lambda x: list(set(x)),
    'score': 'max',
    'dpsi':'max',
    'fisher_pvalue': 'min',
    'junctionLabel': lambda x:   list(set(x)),
    'count_normalize_label': lambda x: list(set(x)),
    }).reset_index()

    integrate_results = result[(result['method'].apply(lambda x: len(x) == 3)) & (result['junctionLabel'].apply(lambda x: len(x) == 1))  & (result['count_normalize_label'].apply(lambda x: len(x) == 1))].reset_index(drop=True)
    integrate_results = integrate_results.explode(['junctionLabel', 'count_normalize_label'])
    overlapped_event = integrate_results['event_label'].unique()

    overlapped_df = df[df['event_label'].isin(overlapped_event)]


    integrate_results['junctionDE_label'] = False
    integrate_results.loc[(integrate_results['fisher_pvalue'] < fisher_pvalue) & 
                        (integrate_results['count_normalize_label'] == True) & 
                        (integrate_results['junctionLabel'] == True),'junctionDE_label'] = True
    integrate_results['psi_label'] = False
    integrate_results.loc[(integrate_results['score'] > score_cutoff) & (integrate_results['dpsi'].abs() > dpsi_cutoff), 'psi_label'] = True
    integrate_results['F1_label'] = integrate_results.apply( lambda x: 'TP' if (x['junctionDE_label'] and x['psi_label'] ) else
                                                                'TN' if (not x['junctionDE_label'] and not x['psi_label']) else
                                                                'FP' if (not x['junctionDE_label'] and x['psi_label']) else
                                                                'FN' if (x['junctionDE_label'] and not x['psi_label']) else 
                                                                'NA',
                                                                axis=1)

    rmats = unique_event_output(overlapped_df[overlapped_df['method'] == 'rmats'])
    leafcutter = unique_event_output(overlapped_df[overlapped_df['method'] == 'leafcutter'])
    majiq = unique_event_output(overlapped_df[overlapped_df['method'] == 'majiq'])

    rmats['junctionDE_label'] = False
    rmats.loc[(rmats['fisher_pvalue'] < fisher_pvalue) & 
                        (rmats['count_normalize_label'] == True) & 
                        (rmats['junctionLabel'] == True),'junctionDE_label'] = True
    rmats['psi_label'] = False
    rmats.loc[(rmats['score'] > score_cutoff) & (rmats['dpsi'].abs() > dpsi_cutoff), 'psi_label'] = True
    rmats['F1_label'] = rmats.apply( lambda x: 'TP' if (x['junctionDE_label'] and x['psi_label'] ) else
                                                                'TN' if (not x['junctionDE_label'] and not x['psi_label']) else
                                                                'FP' if (not x['junctionDE_label'] and x['psi_label']) else
                                                                'FN' if (x['junctionDE_label'] and not x['psi_label']) else 
                                                                'NA',
                                                                axis=1)
    
    leafcutter['junctionDE_label'] = False
    leafcutter.loc[(leafcutter['fisher_pvalue'] < fisher_pvalue) & 
                        (leafcutter['count_normalize_label'] == True) & 
                        (leafcutter['junctionLabel'] == True),'junctionDE_label'] = True
    leafcutter['psi_label'] = False
    leafcutter.loc[(leafcutter['score'] > score_cutoff) & (leafcutter['dpsi'].abs() > dpsi_cutoff), 'psi_label'] = True
    leafcutter['F1_label'] = leafcutter.apply( lambda x: 'TP' if (x['junctionDE_label'] and x['psi_label'] ) else
                                                                'TN' if (not x['junctionDE_label'] and not x['psi_label']) else
                                                                'FP' if (not x['junctionDE_label'] and x['psi_label']) else
                                                                'FN' if (x['junctionDE_label'] and not x['psi_label']) else 
                                                                'NA',
                                                                axis=1)
    
    majiq['junctionDE_label'] = False
    majiq.loc[(leafcutter['fisher_pvalue'] < fisher_pvalue) & 
                        (majiq['count_normalize_label'] == True) & 
                        (majiq['junctionLabel'] == True),'junctionDE_label'] = True
    majiq['psi_label'] = False
    majiq.loc[(majiq['score'] > score_cutoff) & (majiq['dpsi'].abs() > dpsi_cutoff), 'psi_label'] = True
    majiq['F1_label'] = majiq.apply( lambda x: 'TP' if (x['junctionDE_label'] and x['psi_label'] ) else
                                                                'TN' if (not x['junctionDE_label'] and not x['psi_label']) else
                                                                'FP' if (not x['junctionDE_label'] and x['psi_label']) else
                                                                'FN' if (x['junctionDE_label'] and not x['psi_label']) else 
                                                                'NA',
                                                                axis=1)
    
    rmats_F1_metric = pd.DataFrame(rmats['F1_label'].value_counts()).rename(columns={'F1_label': 'rmats'})
    leafcutter_F1_metric = pd.DataFrame(leafcutter['F1_label'].value_counts()).rename(columns={'F1_label': 'leafcutter'})
    majiq_F1_metric = pd.DataFrame(majiq['F1_label'].value_counts()).rename(columns={'F1_label': 'majiq'})
    intergrate_F1_metric = pd.DataFrame(integrate_results['F1_label'].value_counts()).rename(columns={'F1_label': 'SpliceHarmonization'})

    metrics = pd.concat([rmats_F1_metric, leafcutter_F1_metric, majiq_F1_metric,intergrate_F1_metric], axis=1).T
    metrics = metrics.fillna(0)
    metrics['total'] = metrics[['TP', 'TN', 'FP', 'FN']].sum(axis=1)


    metrics['acc'] = metrics.apply(lambda x: (x['TP'] + x['TN']) / (x['TP'] + x['TN'] + x['FP'] + x['FN']) if (x['TP'] + x['TN'] + x['FP'] + x['FN']) != 0 else 0, axis=1)
    metrics['precision'] = metrics.apply(lambda x: x['TP'] / (x['TP'] + x['FP']) if (x['TP'] + x['FP']) != 0 else 0, axis=1)
    metrics['recall'] = metrics.apply(lambda x: x['TP'] / (x['TP'] + x['FN']) if (x['TP'] + x['FN']) != 0 else 0, axis=1)
    metrics['F1'] = metrics.apply(lambda x: 2 * x['precision'] * x['recall'] / (x['precision'] + x['recall']) if (x['precision'] + x['recall']) != 0 else 0, axis=1)
    
    return metrics


def majiq_score_list(row, score_cutoff):
    values = [row['score_cutoff0.1'], row['score_cutoff0.2'], row['score_cutoff0.3'], row['score_cutoff0.4']]
    
    if any(v > score_cutoff for v in values):
        min_value = min(v for v in values if v > score_cutoff)
        index = values.index(min_value)
    else:
        min_value = row['score_cutoff0.1']
        index = 0
    
    return pd.Series([values, min_value, index], index=['score_list', 'score', 'score_index'])

def unique_flattened_list(series):
    # First, evaluate each string in the series to convert it to a list
    evaluated_series = map(ast.literal_eval, series)
    # Flatten the list of lists
    flattened_list = [item for sublist in evaluated_series for item in sublist]
    # Convert the flattened list to a set to get unique items
    unique_set = set(flattened_list)
    return unique_set


def unique_event_df(df):
    res_ = df.groupby(['event_label', 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'annotation': 'first',
            'splice': 'first',
        }
        ).reset_index()

    return res_

def unique_event_df_junction_wopsi(df):
    res_ = df.groupby(['event_label', 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'annotation': 'first',
            'splice': 'first',
            'junctionLabel': lambda x: any(x),
            'count_normalize_label': lambda x: any(x)
        }
        ).reset_index()

    return res_

def unique_event_df_junction(df):
    res_ = df.groupby(['event_label', 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'PSI_1': 'mean',
            'PSI_2': 'mean',
            'annotation': 'first',
            'splice': 'first',
            'junctionLabel': lambda x: any(x),
            'count_normalize_label': lambda x: any(x)
        }
        ).reset_index()

    return res_

def unique_event_df_w_seqs(df):
    res_ = df.groupby(['event_label', 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'PSI_1': 'mean',
            'PSI_2': 'mean',
            'annotation': 'first',
            'splice': 'first',
            'seqs': 'first'
        }
        ).reset_index()

    return res_

def unique_event_df_w_seqs_str(df, str_):
    res_ = df.groupby([str_, 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'PSI_1': 'mean',
            'PSI_2': 'mean',
            'annotation': 'first',
            'splice': 'first',
            'seqs': 'first'
        }
        ).reset_index()

    return res_


def unique_event_df_w_seqs_junction_str(df, str_):
    res_ = df.groupby([str_, 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'PSI_1': 'mean',
            'PSI_2': 'mean',
            'annotation': 'first',
            'splice': 'first',
            'seqs': 'first', 
            'ref_junction_mean': 'mean',
            'alt_junction_mean': 'mean'
        }
        ).reset_index()

    return res_


def unique_event_df_4SpliceAI(df, str_):
    res_ = df.groupby([str_, 'gene_name']).agg(
        {
            'method': lambda x: list(set(x)), 
            'score': 'max',
            'dpsi': 'max',
            'PSI_1': 'mean',
            'PSI_2': 'mean',
            'annotation': 'first',
            'splice': 'first',
            'seqs': 'first', 
            'chr': 'first',
            'strand': 'first'
        }
        ).reset_index()

    return res_

def combine_comparison(data_dict, comparisons):
    subdict = {k: data_dict[k] for k in comparisons}
    df_all = pd.concat(subdict)
    
    return df_all