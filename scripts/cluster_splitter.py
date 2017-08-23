# Script for splitting .rcm file. Creates new .rcm with splitted clusters.
import argparse
import pickle
from sklearn.preprocessing import scale
import func_tools
reload(func_tools)
from func_tools import First_lvl_stacking
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-s', help='source FASTA file with Illumina reads')
parser.add_argument('-r', help='RCM output of IgRepertoireConstructor')
parser.add_argument('-o', help='Output directory')
parser.add_argument('-f', help='Name of the output file')
args = parser.parse_args()

log_file = open(args.o + 'log.txt', 'w')

data = func_tools.build_df(input_reads = args.s, fa_reference = None, rcm_file = args.r)

data = data[['size', 'value1', 'value2', 'value3']]

data['second_vote_abs1'] = data['value1'] * data['size']
data['second_vote_abs2'] = data['value2'] * data['size']
data['second_vote_abs3'] = data['value3'] * data['size']

data = data[['second_vote_abs1', 'second_vote_abs2', 'second_vote_abs3',  'size']]
print 'All clusters in dataset: ', len(data)
log_file.write('All clusters: ' + str(len(data)) + '\n')
data = data[data['size'] > 10]
print 'Clusters with sizes more than 10: ', len(data)
log_file.write('Clusters with sizes more than 10: ' + str(len(data)) + '\n')

target_index = data.index

pp = pickle.load(open('../data/models/pipe_model', 'rb'))

data = pp.steps[0][1].transform(data)
data = pd.DataFrame(data)
data.columns = ['0','1','2','3'] 
logreg = pd.Series(pp.steps[1][1].models['first_lvl'][0].predict(data), index=target_index)
xgb  = pd.Series(pp.steps[1][1].models['first_lvl'][1].predict(data), index=target_index)

ans = pp.steps[2][1].predict(pd.concat([logreg, xgb], axis=1))
#ans = xgb

ans = pd.Series(ans, index = target_index)
target_clusters = ans[ans == 1]

print 'Find ' + str(len(target_clusters)) + ' clusters for splitting in ' + str(len(ans)) + ' clusters'
log_file.write('Find ' + str(len(target_clusters)) + ' clusters for splitting in ' + str(len(ans)) + ' clusters' + '\n')

id_dict = func_tools.id_to_read(args.s)
igrec_rcm = func_tools.read_rcm(args.r)

clusters = func_tools.construct_clusters(igrec_rcm, id_dict)
print "Clusters in IgReC output: ", len(clusters), '\n'
log_file.write('Clusters in IgReC output: ' + str(len(clusters)) + '\n' )

clusters = func_tools.clusters_splitting(clusters, list(target_clusters.index), threshold=5)
print "Clusters in splitted IgReC output: ", len(clusters), '\n'
log_file.write('Clusters in splitted IgReC output: ' + str(len(clusters)) + '\n' )

print 'New ', args.f, '.rcm file consists of ' + str(len(clusters))
log_file.write('New ' + args.f + '.rcm file consists of ' + str(len(clusters)))

func_tools.to_rcm(args.o + args.f + '.rcm', clusters, id_dict)
func_tools.to_fa(args.o + args.f + '.fa', clusters)
log_file.close()
