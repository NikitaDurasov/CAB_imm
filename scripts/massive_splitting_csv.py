import func_tools
import os
import re
import pandas as pd
import pickle

def run_csv_splitter(data, input_reads, igrec_rcm_file, output, filename):
      data['second_vote_abs1'] = data['value1'] * data['size']
      data['second_vote_abs2'] = data['value2'] * data['size']
      data['second_vote_abs3'] = data['value3'] * data['size']

      data = data[['second_vote_abs1', 'second_vote_abs2', 'second_vote_abs3',  'size']]
      print 'All clusters: ', len(data)
      data = data[data['size'] > 10]
      print 'Clusters with sizes more than 10: ', len(data)
      log.write('Clusters with sizes more than 10: ' + str(len(data)) + '\n')
      target_index = data.index

      pp = pickle.load(open('/home/ndurasov/ig_cluster_splitter/data/models/pipe_model', 'rb'))

      data = pp.steps[0][1].transform(data)
      data = pd.DataFrame(data)
      data.columns = ['0','1','2','3'] 
      logreg = pd.Series(pp.steps[1][1].models['first_lvl'][0].predict(data), index=target_index)
      xgb  = pd.Series(pp.steps[1][1].models['first_lvl'][1].predict(data), index=target_index)

      ans = pp.steps[2][1].predict(pd.concat([logreg, xgb], axis=1))

      ans = pd.Series(ans, index = target_index)
      target_clusters = ans[ans == 1]

      print 'Find ' + str(len(target_clusters)) + ' clusters for splitting in ' + str(len(ans)) + ' clusters'
      log.write('Find ' + str(len(target_clusters)) +
              ' clusters for splitting in ' + str(len(ans)) + ' clusters\n')
      id_dict = func_tools.id_to_read(input_reads)

      igrec_rcm = func_tools.read_rcm(igrec_rcm_file)

      clusters = func_tools.construct_clusters(igrec_rcm, id_dict)
      clusters = func_tools.clusters_splitting(clusters, list(target_clusters.index))

      print 'New .rcm file consists of ' + str(len(clusters))
      log.write('New .rcm file consists of ' + str(len(clusters)) + '\n\n')

      func_tools.to_rcm(output + filename + '.rcm', clusters, id_dict)
      func_tools.to_fa(output + filename + '.fa', clusters)

base = '/Bmo/orange_nikita/'
directories =  list(os.walk(base))[0][1]

regex = re.compile('.*_compilation')
compilations = filter(regex.match, directories)

for compilation in compilations:
    if not os.path.exists(base + compilation + '/splitted'):
            os.makedirs(base + compilation + '/splitted')
    regex_csv = re.compile('.*\.csv')
    igrec_csv = filter(regex_csv.match, os.listdir(base + compilation))
    if 'compiled_dataset.csv' in igrec_csv:
        igrec_csv.remove('compiled_dataset.csv')
    for dataset in igrec_csv:
        log = open(base + compilation + '/splitted/log.txt', 'a')
        print 'Splitting ' + dataset
        igrec_output = dataset[:-12]
        df = pd.read_csv(base + compilation + '/' + dataset) 
        log.write('Splitting ' + dataset + '\n')
        log.write('Lenght: ' + str(len(df)) + '\n')
        input_reads = base + compilation + '/' + igrec_output + '/input3.fa'
        igrec_rcm = base + compilation + '/' + igrec_output  + '/final_repertoire.rcm'
        output = base + compilation + '/splitted/'
        filename = igrec_output + '_splitted'
        run_csv_splitter(df, input_reads, igrec_rcm, output, filename)
        log.close()


