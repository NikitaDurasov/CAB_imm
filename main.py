from func_tools import *
import pickle
from Bio import SeqIO

with open('final_clusters.pkl', 'rb') as inputio:
    final_clusters = pickle.load(inputio)
with open('final_rep.pkl', 'rb') as inputio:
    final_rep = pickle.load(inputio)
with open('final_res.pkl', 'rb') as inputio:
    final_res = pickle.load(inputio)

# path = "/Users/Macbook/GitHub/ig_cluster_splitter/test/"
# build_df(input_reads = path + 'input_reads_test.fa',
#         fa_reference = path + 'final_repertoire.fa',
#          rcm_file = path + 'final_repertoire.rcm')

build_df_test(final_clusters, final_rep, final_res)

