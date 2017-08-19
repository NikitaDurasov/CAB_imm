# Collect info into compilation directory  
import re
import subprocess as sb
import os
from shutil import copy

regex_str = raw_input('Enter regex for directories: ')

repertoire_path = '/Bmo/ashlemov/ig_repertoire_constructor_jul5/src/extra/ig_quast_tool'

all_dirs = os.walk('..').next()[1]
regex = re.compile(str(regex_str))
regexed_dirs = filter(regex.match, all_dirs)
curr_dir = os.getcwd().split('/')[-1]

if curr_dir in regexed_dirs:
    regexed_dirs.remove(os.getcwd().split('/')[-1])

for dir in regexed_dirs:
	res_fa = '../' + dir + '/' + 'final_repertoire.fa' 
	res_rcm = '../' + dir + '/' + 'final_repertoire.rcm'
	input_fa = repertoire_path + '/' + dir[:-13]  + '/' + 'input3.fa.gz'
	ref_fa = repertoire_path + '/' + dir[:-13]  + '/' + 'repertoire3.fa.gz'
	ref_rcm = repertoire_path + '/' + dir[:-13]  + '/' + 'repertoire3.rcm'
	
	dir_short = dir[:-13]

	destination_directory = os.getcwd() + '/' +  dir_short
	
	if not os.path.exists(destination_directory):
		print 'Create directory: ' + destination_directory 
		os.mkdir(destination_directory)

	if not os.path.exists(destination_directory + 'final_repertoire.fa'):
		print 'Copy ' + res_fa + ' to ' + destination_directory
		copy(res_fa, dir_short)
	else:
		os.remove(destination_directory + 'final_repertoire.fa')

	if not os.path.exists(destination_directory + 'final_repertoire.rcm'):
		print 'Copy ' + res_rcm + ' to ' + destination_directory
		copy(res_rcm, dir_short)
	else:
		os.remove(destination_directory + 'final_repertoire.rcm')

	if not os.path.exists(destination_directory + 'input3.fa.gz'):
		print 'Copy ' + input_fa + ' to ' + destination_directory
		copy(input_fa, dir_short)
	else:
		os.remove(destination_directory + 'input3.fa.gz')

	if not os.path.exists(destination_directory + 'repertoire3.fa.gz'):
		print 'Copy ' + ref_fa + ' to ' + destination_directory
		copy(ref_fa, dir_short)
	else:
		os.remove(destination_directory + 'repertoire3.fa.gz')

	if not os.path.exists(destination_directory + 'repertoire.rcm'):
		print 'Copy ' + res_rcm + ' to ' + destination_directory
		copy(ref_rcm, dir_short)
	else: 
		os.remove(destination_directory + 'repertoire.rcm')




