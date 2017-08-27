import func_tools
reload(func_tools)
import os
import re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def draw_sizes_hist(filename, a, b, dataset):
    fig = plt.figure(figsize=(30, 10))
    plt.grid(alpha=0.3)
    plt.bar(range(5, 50), a[5:50], alpha=0.5, color='y', label='IgReC output')
    plt.bar(range(5, 50), b[5:50], alpha=0.6, color='r', label='Stacking model')
    plt.xticks(range(5, 50, 5), range(5, 50, 5), fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim((4, 36))
    plt.title(dataset, fontsize=20)
    plt.legend(prop={'size':30})

    fig.savefig(filename, format='pdf')

base = '/Bmo/orange_nikita/general/august24/'
files = os.listdir(base)
regex = re.compile('.*\.rcm')
files = filter(regex.match, files)

dataset_directory = {}
dataset_directory['AGE'] = '/Bmo/orange_nikita/AGE_compilation/'
dataset_directory['FLU'] = '/Bmo/orange_nikita/FLU_compilation/'
dataset_directory['GMC'] = '/Bmo/orange_nikita/GMC_compilation/'
dataset_directory['IDO'] = '/Bmo/orange_nikita/IDO_compilation/'

for rcm_file in sorted(files):

    curr_dir = dataset_directory[rcm_file[:3]]
    curr_dataset = rcm_file.split('splitted')[0][:-1]
    input_reads = curr_dir + curr_dataset + '/input3.fa'
    ref_rcm_file = curr_dir + curr_dataset + '/repertoire3.rcm'
    cons_rcm_file_igrec = curr_dir + curr_dataset + '/final_repertoire.rcm'
    cons_rcm_file_splitted = base + rcm_file

    igrec_dict = func_tools.unrecognized_clusters_sizes(input_reads, ref_rcm_file, cons_rcm_file_igrec)
    splitted_dict = func_tools.unrecognized_clusters_sizes(input_reads, ref_rcm_file, cons_rcm_file_splitted)

    filename_perc_sizes = base + 'graphs/' + curr_dataset  + '_perc_sizes.pdf'
    print 'Drawing ' + curr_dataset + '_perc_sizes.pdf'
    draw_sizes_hist(filename_perc_sizes, igrec_dict['perc_sizes_hist'], splitted_dict['perc_sizes_hist'], curr_dataset)
    print 'File saved to ' + filename_perc_sizes
