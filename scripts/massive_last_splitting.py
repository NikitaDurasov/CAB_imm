import func_tools
import os
import re

base = '/Bmo/orange_nikita/'
directories = list(os.walk(base))[0][1]
output = '/Bmo/orange_nikita/general/august24/'

regex = re.compile('.*_compilation')
compilations = filter(regex.match, directories)

for compilation in compilations:
    curr_compilation_directory = base + compilation
    datasets = next(os.walk(curr_compilation_directory))[1]

    if 'splitted' in datasets:
        datasets.remove('splitted')

    if 'sizes_hists' in datasets:
        datasets.remove('sizes_hists')

    for dataset in datasets:

        try:
            input_reads = curr_compilation_directory + '/' + dataset + '/input3.fa'
            ref_rcm = curr_compilation_directory + '/' + dataset + '/final_repertoire.rcm'
            output_dir = output
            filename = dataset + '_splitted_5_thesh_1_model'

            print 'Splitting ' + dataset + ' into directory ' + output_dir + ' with name ' + filename + '\n' 

            command = 'python cluster_splitter.py'
            command += ' -s ' + input_reads
            command += ' -r ' + ref_rcm
            command += ' -o ' + output_dir
            command += ' -f ' + filename

            os.system(command)

        except IOError as e:
            log_file = open(output + 'err_log.txt', 'w')
            print >> log_file, str(e)
            print str(e)





