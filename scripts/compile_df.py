# Scipt for completed datasets compilation in base directory  
#
import pandas as pd
import os
import re

def compile_df(base):

    df_arr = []

    files = os.listdir(base)
    regex = re.compile('.+_2\.csv')
    csv_files = filter(regex.match, files)

    if 'compiled_dataset.csv' in csv_files:
        csv_files.remove('compiled_dataset.csv')

    csv_files = map(lambda x: base + x, csv_files)

    print 'files for compilation: '
    for i in csv_files:
        print i

    for csv_file in csv_files:
            temp_df = pd.read_csv(csv_file)
            temp_df['dataset'] = csv_file[:-4]
            df_arr.append(temp_df)
    full_df = pd.concat(df_arr, ignore_index = True)
    full_df.to_csv(base + 'compiled_dataset.csv', index = False)


if __name__ == '__main__':

    base = input("base for script: ")

    df_arr = []

    files = os.listdir(base)
    regex = re.compile('.+\.csv')
    csv_files = filter(regex.match, files)

    if 'compiled_dataset.csv' in csv_files:
        csv_files.remove('compiled_dataset.csv')

    print 'files for compilation: '
    for i in csv_files:
        print i

    for csv_file in csv_files:
            temp_df = pd.read_csv(csv_file)
            temp_df['dataset'] = csv_file[:-4]
            df_arr.append(temp_df)
    full_df = pd.concat(df_arr, ignore_index = True)
    full_df.to_csv('compiled_dataset.csv', index = False)
