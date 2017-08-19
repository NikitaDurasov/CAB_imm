#!/usr/bin/python2.7
import subprocess as sb
import os
import re

dir_path = '/Bmo/ashlemov/ig_repertoire_constructor/src/extra/ig_quast_tool'
igrec_tool_path = '/home/ndurasov/ig_repertoire_constructor/ig_repertoire_constructor'
output_dir_path = '/Bmo/orange_nikita'

files_arr = os.listdir(dir_path)
regex = re.compile(r'IDO_([\d]*)_IG+$')
dirs = filter(regex.match, files_arr)

output_dirs = []
command = ''

for dir in dirs:
	output_dirs.append(dir  + '_igrec_output')

for output_dir, dir  in zip(output_dirs, dirs):
        print "Running igrec on: ", output_dir
        print '\n'
	o  =  output_dir_path + '/' +  output_dir
	s = dir_path + '/' + dir + '/' + 'input3.fa.gz'
	args = '-o ' + o + ' -s ' + s + ' --no-alignment'
	command = igrec_tool_path + '/' + 'igrec.py ' + args
	os.system('mkdir  ' + o)
	os.system(command)

print 'OK'

