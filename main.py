from func_tools import *
import pickle
from Bio import SeqIO

#with open('final_clusters.pkl', 'rb') as inputio:
#    final_clusters = pickle.load(inputio)
#with open('final_rep.pkl', 'rb') as inputio:
#    final_rep = pickle.load(inputio)
#with open('final_res.pkl', 'rb') as inputio:
#    final_res = pickle.load(inputio)

with open('ans_dict.pkl', 'rb') as inputio:
    ans_dict = pickle.load(inputio)


build_ans_df(ans_dict)

