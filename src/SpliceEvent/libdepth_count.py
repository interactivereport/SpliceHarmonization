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
from multiprocessing import Pool

# Set the environment variables for the CA bundles
os.environ['REQUESTS_CA_BUNDLE'] = '/etc/pki/ca-trust/extracted/openssl/ca-bundle.trust.crt'
os.environ['CURL_CA_BUNDLE'] = '/etc/pki/ca-trust/extracted/openssl/ca-bundle.trust.crt'

parser = argparse.ArgumentParser()
parser.add_argument('-outdir', dest="outdir", help="file path")
parser.add_argument('-Altbam', dest='Altbam', type=str, help="file path list")
parser.add_argument('-Refbam', dest='Refbam', type=str, help="file path list")
args = parser.parse_args()

print(args.Altbam)
print(args.Refbam)


# def get_counts(bam_file_paths, threads=8):
#                 with Pool(threads) as pool:
#                                 results = pool.map(bam_count, bam_file_paths)
#                                 # assumes each file path is unique in the list
#                                 results = {result[1]: result[0] for result in results}
#                 return results  


# def bam_count(bam_file_path):
#                 with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
#                                 # return a tuple for traceability in async processes
#                                 return (round(bam_file.count()/1e6), bam_file_path)

try:
    alt_bam_files = ast.literal_eval(args.Altbam)
    ref_bam_files = ast.literal_eval(args.Refbam)

    alt_depth = process_bam_count(alt_bam_files)
    ref_depth = process_bam_count(ref_bam_files)
    # alt_depth = get_counts(alt_bam_files)
    # ref_depth = get_counts(ref_bam_files)

    alt_depth_dict = dict(zip(alt_bam_files, alt_depth))
    ref_depth_dict = dict(zip(ref_bam_files, ref_depth))

except ValueError:
    print("Invalid format. Please provide a valid list of file paths.")
    print(alt_bam_files)
except SyntaxError:
    print("Syntax error. Please provide a proper list format.")


depth_dict = {'alt_depth':alt_depth_dict, 'ref_depth': ref_depth_dict}
print(depth_dict)


file_name = args.outdir + '/libdepth_counts.json'

# Writing JSON data
with open(file_name, 'w') as f:
    json.dump(depth_dict, f, indent=4)