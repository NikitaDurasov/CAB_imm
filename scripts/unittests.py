import unittest
import func_tools
from collections import Counter
reload(func_tools)

test_base = '/home/ndurasov/ig_cluster_splitter/scripts/test_data/'

test_rcm_file = 'test_rcm.rcm'
test_rcm = {'id_number_4': '3',
             'id_number_3': '3',
             'id_number_2': '2',
             'id_number_1': '1'}

test_fa_file = 'test_fa.fa'
test_repertoire = ['CAGCTGC', 'TCCTGCAATG', 'CTGCAT', 'CTGCAA']

test_input_reads_file = 'test_input.fa'
test_input = {'id_number_4': 'CTGCAA',
              'id_number_3': 'CTGCAT',
              'id_number_2': 'TCCTGCAATG',
              'id_number_1': 'CAGCTGC'}

test_clusters = {'1': ['CAGCTGC'], '3': ['CTGCAA', 'CTGCAT'], '2': ['TCCTGCAATG']} 

test_variety = {5: {'A': 1, 'T': 1}}

test_second_vote = {5: 0.5}

class TestReadingFunctions(unittest.TestCase):

    def test_read_rcm(self):
        self.assertEqual(func_tools.read_rcm(test_base + test_rcm_file), test_rcm)

    def test_read_repertoire(self):
        self.assertEqual(func_tools.read_repertoire(test_base + test_fa_file), test_repertoire)

    def test_id_to_read(self):
        self.assertEqual(func_tools.id_to_read(test_base + test_input_reads_file), test_input)

    def test_construct_clusters(self):
       rcm = func_tools.read_rcm(test_base + test_rcm_file)
       input_reads = func_tools.id_to_read(test_base + test_input_reads_file)
       self.assertEqual(func_tools.construct_clusters(rcm, input_reads), test_clusters)


class TestClustersFunctions(unittest.TestCase):

    @classmethod
    def setUpClass(self):
       self.rcm = func_tools.read_rcm(test_base + test_rcm_file)
       self.input_reads = func_tools.id_to_read(test_base + test_input_reads_file)
       self.clusters = func_tools.construct_clusters(self.rcm, self.input_reads)

    def test_cluster_variety(self):
       self.assertEqual(func_tools.cluster_variety(self.clusters['3']), test_variety) 
       self.assertEqual(func_tools.cluster_variety(self.clusters['1']), {})

    def test_second_vote(self):
       self.assertEqual(func_tools.second_vote(self.clusters['1']), {})
       self.assertEqual(func_tools.second_vote(self.clusters['3']), test_second_vote)

    def test_major_vote(self):
        self.assertEqual(func_tools.major_vote(self.clusters['3']), 'CTGCAT')
        self.assertEqual(func_tools.major_vote(self.clusters['1']), 'CAGCTGC')

    def test_second_vote_letter(self):
        self.assertEqual(func_tools.second_vote_letter(self.clusters['3'], 5), 'A')
        self.assertEqual(func_tools.second_vote_letter(self.clusters['3'], 6), ' ')
        self.assertEqual(func_tools.second_vote_letter(self.clusters['3'], 4), ' ')

    def test_clusters_size_dict(self):
         self.assertEqual(func_tools.clusters_size_dict(self.clusters), {'1': 1, '3': 2, '2': 1})

    def test_split_cluster(self):
        self.assertEqual(func_tools.split_cluster(self.clusters['3'], 5), ([], []))
        self.assertEqual(func_tools.split_cluster(self.clusters['3'], 5, threshold=0),(['CTGCAA'], ['CTGCAT']))
        self.assertEqual(func_tools.split_cluster(self.clusters['1'], 2), ([],[]))

    def test_find_max_dict(self):
        self.assertEqual(func_tools.find_max_dict({'one':1, 'two':2, 'three':3, 'four':4}), 'four')

    def test_check_size(self):
        self.assertEqual(func_tools.check_size(self.clusters['1'],self.clusters['2']), ([],[]))
        self.assertEqual(func_tools.check_size(self.clusters['1'],self.clusters['2'], threshold=0), 
                         (self.clusters['1'], self.clusters['2']))

    def test_split_by_2nd_vote(self):
        self.assertEqual(func_tools.split_by_2nd_vote(self.clusters['3'], threshold=0), (['CTGCAA'], ['CTGCAT']))
        self.assertEqual(func_tools.split_by_2nd_vote(self.clusters['3']), ([],[]) )

    def test_clusters2rep(self):
        self.assertEqual(func_tools.clusters2rep(self.clusters), {'1': 'CAGCTGC', '3': 'CTGCAT', '2': 'TCCTGCAATG'})
        self.assertEqual(func_tools.clusters2rep({}), {})

    def test_clusters_filtering(self):
        self.assertEqual(func_tools.clusters_filtering(self.clusters), {})
        self.assertEqual(func_tools.clusters_filtering(self.clusters, threshold=0), self.clusters)
        temp_clusters = {'1':['AC','AC','AC','AC','AC','AC']}
        self.assertEqual(func_tools.clusters_filtering(temp_clusters), temp_clusters)

    def test_second_votes(self):
        self.assertEqual(func_tools.second_votes({'1': ['AAA', 'AAA', 'AA']}), {'1': {}})
        self.assertEqual(func_tools.second_votes({'1': ['AAA', 'AAA', 'AAC']}), {'1': {2: 1.0 / 3}})
        self.assertEqual(func_tools.second_votes({}), {})



if __name__ == '__main__':
    unittest.main()
