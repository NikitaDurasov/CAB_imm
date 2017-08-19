import os 
base = '/Bmo/orange_nikita/AGE_compilation'

all_dirs = sorted(os.walk('.').next()[1])

for dir in all_dirs:
	print 'Unzip files in ' + dir 
	os.system('gunzip -k ' + dir + '/' + '*.gz')
