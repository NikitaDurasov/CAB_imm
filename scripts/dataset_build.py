# script for building dataset from .fa and .rcm files

import func_tools as ft
reload(ft)
import pandas as pd
import re
import subprocess as sb
import os

def dataset_build(base, suffix):
    regex_str = raw_input('Enter regex for directories: ')

    all_dirs = os.walk(base).next()[1]
    regex = re.compile(str(regex_str))
    regexed_dirs = filter(regex.match, all_dirs)

    input_reads = [base + dir + '/input3.fa' for dir in regexed_dirs]
    references = [base + dir + '/repertoire3.fa' for dir in regexed_dirs]
    reference_rcms = [base + dir + '/repertoire3.rcm' for dir in regexed_dirs]
    final_rcms = [base + dir + '/final_repertoire.rcm' for dir in regexed_dirs]

    for input, ref, final, dir, ref_rcm  in zip(input_reads, references, final_rcms, regexed_dirs, reference_rcms):
        print 'Building dataset ' + dir
        temp_df = ft.build_df(input, final, ref_rcm, ref, classification='colormap', threshold=0)
        temp_df['dataset'] = [dir] * len(temp_df)
        temp_df.to_csv(base + '/' + dir + '_dataset' + suffix + '.csv')

if __name__ == '__main__':

    base = input('Enter base direcory: ')
    suffix = input('Enter suffix for files: ')
    regex_str = raw_input('Enter regex for directories: ')

    all_dirs = os.walk(base).next()[1]
    regex = re.compile(str(regex_str))
    regexed_dirs = filter(regex.match, all_dirs)

    input_reads = [base + dir + '/input3.fa' for dir in regexed_dirs]
    references = [base + dir + '/repertoire3.fa' for dir in regexed_dirs]
    final_rcms = [base + dir + '/final_repertoire.rcm' for dir in regexed_dirs]

    for input, ref, final, dir in zip(input_reads, references, final_rcms, regexed_dirs):
        print 'building flu' + dir
        temp_df = ft.build_df(input, ref, final)
        temp_df['dataset'] = [dir] * len(temp_df)
        temp_df.to_csv(base + '/' + dir + '_dataset' + suffix + '.csv')
