import func_tools as ft
import os
import re

def run_cluster_splitter(input_reads, igrec_rcm, output, filename):
    splitter  = '/home/ndurasov/ig_cluster_splitter/scripts/cluster_splitter.py'
    os.system('python ' + splitter + ' -s ' + input_reads +
                ' -r ' + igrec_rcm + ' -o ' + output +  ' -f ' + filename)

base = '/Bmo/orange_nikita/'
directories =  list(os.walk(base))[0][1]

regex = re.compile('.*_compilation')
compilations = filter(regex.match, directories)

for compilation in compilations:
    if not os.path.exists(base + compilation + '/splitted'):
            os.makedirs(base + compilation + '/splitted')
    igrec_outputs =  list(os.walk(base + compilation))[0][1]
    igrec_outputs.remove('splitted')
    for igrec_output in igrec_outputs:
        print 'Splitting ' + igrec_output
        input_reads = base + compilation + '/' + igrec_output + '/input3.fa'
        igrec_rcm = base + compilation + '/' + igrec_output  + '/final_repertoire.rcm'
        output = base + compilation + '/splitted/'
        filename = compilation + '_splitted.rcm'
        run_cluster_splitter(input_reads, igrec_rcm, output, filename)
