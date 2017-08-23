import func_tools
reload(func_tools)
import os

base = '/Bmo/orange_nikita/'
directories = list(ps.walk(base))[0][1]

regex - re.compile('.*_compilation')
compilations = filter(regex.match, directories)

for compilation in compilations:
    curr_dir = base + compilation + '/sizes_hists'

    if not os.path.exists(curr_dir):
        os.makedir(curr_dir)

#FINISH THIS SHIT 



