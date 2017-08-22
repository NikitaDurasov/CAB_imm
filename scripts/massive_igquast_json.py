import func_tools
import os
import re

igquast_loc = '/home/ndurasov/ig_repertoire_constructor/ig_repertoire_constructor/igquast.py'
base = '/Bmo/orange_nikita/'
directories =  list(os.walk(base))[0][1]

regex = re.compile('.*_compilation')
compilations = filter(regex.match, directories)

for compilation in compilations:
    curr_dir = base + compilation + '/splitted'
    regex_fa = re.compile('.*\.fa+$')
    fa_files = filter(regex_fa.match, os.listdir(curr_dir))
    fa_files = sorted(fa_files)
    for fa in fa_files:
        dataset = fa[:-12]
        if not os.path.exists(curr_dir + '/' + dataset + '_igquast/igrec_splitted'):
            os.makedirs(curr_dir + '/' + dataset + '_igquast/igrec_splitted')
        if not os.path.exists(curr_dir + '/' + dataset + '_igquast/igrec_origin'):
            os.makedirs(curr_dir + '/' + dataset + '_igquast/igrec_origin')

        input_reads = base + compilation + '/' + dataset + '/input3.fa'
        reference = base + compilation + '/' +  dataset + '/repertoire3.fa'

        igrec_origin = base + compilation + '/' + dataset + '/final_repertoire.fa'
        splitted = curr_dir + '/' + fa
        output_splitted = curr_dir + '/' + dataset + '_igquast/igrec_splitted'
        output_origin = curr_dir + '/' + dataset + '_igquast/igrec_origin'

        os.system(igquast_loc + ' -s ' + input_reads + ' -r ' +
                    reference + ' -c  ' + splitted + ' -o ' + output_splitted)
        os.system(igquast_loc + ' -s ' + input_reads + ' -r ' +
                    reference + ' -c  ' + igrec_origin + ' -o ' + output_origin)

