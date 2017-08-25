import func_tools
reload(func_tools)
import os
import re
import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#def draw_sizes_hist(filename, arr1, arr2, dataset):
#    fig = plt.figure(figsize=(30, 10))
#    plt.grid(alpha=0.3)
#    plt.bar(range(5,75), arr1[5:75], alpha=0.5, color='y', label='IgReC output')
#    plt.bar(range(5,75), arr2[5:75], alpha=0.6, color='r', label='Stacking model')
#    plt.xticks(range(5, 75, 5), range(5, 75, 5), fontsize=20)
#    plt.yticks(fontsize=20)
#    plt.xlim((4, 35))
#    plt.title(dataset, fontsize=20)
#    plt.legend(prop={'size':30})
#    fig.savefig(filename, format='pdf')

def draw_sizes_hist(filename, a, b, dataset)
    plt.figure(figsize=(30, 10))
    plt.grid(alpha=0.3)
    plt.bar(range(5,200), a[:195], alpha=0.5, color='y', label='IgReC output')
    plt.bar(range(5,200), b[:195], alpha=0.6, color='r', label='Stacking model')
    plt.xticks(range(5, 200, 5), range(5, 200, 5), fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim((4, 30))
    plt.title(dataset, fontsize=20)
    plt.legend(prop={'size':30})

    fig.savefig(filename, format='pdf')

base = '/Bmo/orange_nikita/'
directories = list(os.walk(base))[0][1]

regex = re.compile('.*_compilation')
compilations = filter(regex.match, directories)

for compilation in compilations:
    curr_dir = base + compilation + '/sizes_hists'

    if not os.path.exists(curr_dir):
        os.mkdir(curr_dir)
        print curr_dir, ' created' 

    subdirectories = next(os.walk(base + compilation))[1]

    if 'sizes_hists' in subdirectories:
        subdirectories.remove('sizes_hists')  
    if 'splitted' in subdirectories:
        subdirectories.remove('splitted')  


    for directory in subdirectories:
        input_reads = base + compilation + '/' + directory + '/input3.fa'
        ref_rcm_file = base + compilation + '/' + directory + '/repertoire3.rcm'
        cons_rcm_file_igrec =  base + compilation + '/' + directory + '/final_repertoire.rcm'
        cons_rcm_file_splitted = base + compilation + '/splitted/' + directory + '_splitted.rcm'

        igrec_dict = func_tools.unrecognized_clusters_sizes(input_reads, ref_rcm_file, cons_rcm_file_igrec)
        splitted_dict = func_tools.unrecognized_clusters_sizes(input_reads, ref_rcm_file, cons_rcm_file_splitted)

        filename_perc_sizes = curr_dir + '/' + directory + '_perc_sizes.pdf'
        print 'Drawing ' + directory + '_perc_sizes.pdf'
        draw_sizes_hist(filename_perc_sizes, igrec_dict['perc_sizes_hist'], splitted_dict['perc_sizes_hist'], directory)
        print 'File saved to ' + filename_perc_sizes
