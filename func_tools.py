import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter, defaultdict, OrderedDict
import nltk
from nltk.util import ngrams
import re
import itertools
import operator

def word_grams(words, min_v=1, max_v=4):
    """
    :param words: one long string without spaces
    :param min_v: min_v length of ngram that'll be returned
    :param max_v: max-1 length of ngram that'll be returned
    :return: array of all ngrams from min_v to max_v-1 collected form words
    """
    s = []
    for n in range(min_v, max_v):
        for ngram in ngrams(words, n):
            s.append(' '.join(str(i) for i in ngram))
    return s


def build_letter_hist(rep_ngram, figsize=(5, 5), threshold=0):
    """
    Draw ngrams frequencies histogram
    :param rep: array of ngrams for every element in repertoire. Example: ['AAAC', 'CCC'] -> rep = [['A A', 'A A', 'A C'], ['C C', 'C C']]
    :param figsize: matplotlib figsize arg for plot
    :param threshold: drop ngrams with frequency lower than threshold
    :return: None
    """
    letter_counts = Counter(np.hstack(rep_ngram))
    letter_counts_new = drop_rare(letter_counts, threshold=threshold)
    df = pd.DataFrame.from_dict(letter_counts_new, orient='index')
    df.plot(kind='bar', figsize=figsize)


def merge_subarrays(arr):
    """
    Merge subarrays in one long array. Example: [[1,2,3], [2,3,4]] -> [1,2,3,2,3,4]
    :param arr: array of arrays to merge
    :return: merged array
    """
    return np.hstack(arr)


def drop_rare(counter, threshold=0):
    """
    Drop elements from counter/dict with frequencies lower than threshold
    :param counter:  Counter or dict object
    :param threshold: int value
    :return: changed Counter or dict object
    """
    for k in list(counter):
        if counter[k] < threshold:
            del counter[k]
    return counter


def calculate_prob(words, n=2):
    """
    Function calculates probability evaluations for letters following ngrams.
    :param words: one long string obj
    :param n: length of ngrams for which this function is used
    :return: dictionary of dicts; ngram -> dict of probabilities 
    """
    temp = np.array(list(words))
    result = {}
    for ngram in set(word_grams(words, min_v=n, max_v=n + 1)):
        ngram = ngram.replace(" ", "")
        ngram_positions = [m.start() for m in re.finditer('(?=' + ngram + ')', words)]
        letters_positions = np.array(ngram_positions) + len(ngram)
        letters_position = letters_positions[letters_positions < len(words)]
        counter = dict(Counter(temp[letters_position]))
        total = sum(counter.itervalues(), 0.0)
        counter = {k: 1.0 * v / total for k, v in counter.iteritems()}
        result[ngram] = counter
    return result


def id_to_read(filename):
    """
    Creates dict object which map read id's from .fa file to actual read sequence
    :param filename: filename of .fa file with read's and id's
    :return: dict id -> sequence
    """
    inp_seq = SeqIO.parse(filename, "fasta")
    inp_seq = list(inp_seq)
    read_id = [x.id for x in inp_seq]
    inp_reads = [str(x.seq) for x in inp_seq]
    id_to_read = {k: v for k, v in zip(read_id, inp_reads)}
    return id_to_read


def read_repertoire(filename):
    """
    Read repertoire from .fa into array
    :param filename: filename of .fa file with repertoire
    :return: array of strings
    """
    return map(lambda x: str(x.seq), SeqIO.parse('repertoire.fa', "fasta"))


def read_rcm(filename):
    """
    Read .rcm file into dict obj "read_id" -> "cluster_id
    :param filename:
    :return:
    """
    rcm = open(filename)
    rcm = rcm.readlines()
    rcm = [x[:-1] for x in rcm]
    return dict(map(lambda x: x.split("\t"), rcm))


def construct_clusters(rcm_dict, id_dict):
    """
    Function constructs clusters dict (cluster number) -> (array of reads included into this cluster)
    :param rcm_dict: dict obj with reads ids -> cluster number
    :param id_dict: dict obj with reads ids -> array of reads included into this cluster
    :return: dict obj with cluster number -> array of reads included into this cluster
    """
    clusters = defaultdict(lambda: [])
    for value, key in rcm_dict.items():
        clusters[key].append(id_dict[value])
    return clusters


def max_character(arr):
    """
    Find maximum element in the string arr (order is desc: a<z)
    :param arr: string or list of characters
    :return: maximum element in the string or list of characters
    """
    max_ch = 'A'
    for letter in arr:
        if max_ch < letter:
            max_ch = letter
    return max_ch


def min_character(arr):
    """
    Find minimum element in the string arr (order is desc: a<z)
    :param arr: string or list of characters
    :return: minimum element in the string or list of characters
    """
    min_ch = 'Z'
    for letter in arr:
        if min_ch > letter:
            min_ch = letter
    return min_ch


def cluster_variety(cluster):
    """
    Constructs default dict obj, that consists of elements (position in cluster) -> Counter obj, where Counter count
    frequencies of every letter appeared in particulat position
    :param cluster: list of reads
    :return: dict odj with (position in cluster) -> Counter obj
    """
    letter_matrix = []
    res = defaultdict(lambda: 0)
    for read in cluster:
        letter_matrix.append(list(read))
    letter_matrix = list(itertools.izip_longest(*letter_matrix))
    for i, row in enumerate(letter_matrix):
        row = [x for x in row if x is not None]
        if max_character(row) != min_character(row):
            res[i] = Counter(row)
    return res


def second_vote(cluster):
    """
    For every position in cluster function calculates second voter value - frequency of second the most frequent letter
    in position
    :param cluster: list of reads
    :return: dict odj with (position in cluster) -> second vote value
    """
    res = defaultdict(lambda: 0)
    for key, value in cluster_variety(cluster).items():
        temp = sorted(value.values())
        res[key] = 1.0 * temp[-2] / len(cluster)
    return res


def major_vote(cluster):
    """
    Calculates word with the most frequent letter in every position
    :param cluster: list of reads
    :return: string
    """
    ans = list(max(cluster, key=len))
    res = {}
    for key, value in cluster_variety(cluster).items():
        temp = sorted(value.items(), key=lambda x: x[1])
        res[key] = temp[-1][0]
    for key in res:
        ans[key] = res[key]
    return ''.join(ans)


def second_votes(clusters):
    """
    For every cluster in clusters and for every position in cluster function calculates second voter value -
    frequency of second the most frequent letter in position
    :param clusters: dict obj with (cluster number) -> (list of reads)
    :return: dict obj with (cluster number) -> (dict obj with (position in cluster) -> second vote value in this position)
     """
    res = {}
    for key in clusters:
        res[key] = second_vote(clusters[key])
    return res


def clusters_size_dict(clusters):
    """
    Constructs dict obj consists of list lengths for every cluster in clusters
    :param clusters: dict obj with (cluster number) -> (list of reads)
    :return: dict obj with (cluster number) -> (list length)
    """
    clusters_lenghts = {}
    for key in clusters:
        clusters_lenghts[key] = len(clusters[key])
    return clusters_lenghts


def top_massive_clusters(clusters, res_dict, n=10):
    """
    Draw second votes graph of first n the most large clusters in clusters
    :param clusters: dict obj with (cluster number) -> (list of reads)
    :param res_dict: dict obj returned by second_votes function
    :param n: number of drawn pictures
    :return: list of list consists of second vote values for every position in n the most large clusters
    """
    clusters_lenghts = clusters_size_dict(clusters)
    top_n = [x[0] for x in sorted(clusters_lenghts.items(), key=lambda x: x[1])[-n:]]
    res = []
    fig, axes = plt.subplots(ncols=5, nrows=(n - 1) / 5 + 1, figsize=(20, 10))
    for j, key in enumerate(top_n):
        res1 = [0] * 500
        for i in range(500):
            res1[i] += res_dict[str(key)][i] / n
        res.append(res1[:])
        axes[j / 5, j % 5].scatter(range(500), res1)
        axes[j / 5, j % 5].set_ylim((0, 0.008))

    return res

#--------------------------------------------------------------------------------------------------------------------

def precision_sensetivity_F1(constructed, reference):
    ref = set(reference)
    con = set(constructed)

    ref_con_intersection = ref.intersection(con)

    precision = 1.0 * len(ref_con_intersection) / len(con)
    sensitivity = 1.0 * len(ref_con_intersection) / len(ref)

    return 2.0 * precision * sensitivity / (precision + sensitivity)


def true_false_count(constructed, reference):
    true = false = 0
    ref = set(reference)
    con = set(constructed)
    for item in con:
        if item in ref:
            true += 1
        else:
            false += 1
    return (true, false)


def split_cluster(cluster, position):
    """
    Try some other variants
    :param cluster:
    :param position:
    :return:
    """

    position_variation = cluster_variety(cluster)[position]
    if not position_variation:
        return [cluster]
    letter_for_split = find_max_dict(position_variation)
    first_splitted_part = []
    second_splitted_part = []

    for item in cluster:
        if item[position] == letter_for_split:
            first_splitted_part.append(item)
        else:
            second_splitted_part.append(item)

    return (first_splitted_part, second_splitted_part)


def find_max_dict(dictionary):
    return max(dictionary.iteritems(), key=operator.itemgetter(1))[0]

def split_by_2nd_vote(cluster):
    max_2nd_vote_pos = find_max_dict(second_vote(cluster))
    return split_cluster(cluster, max_2nd_vote_pos)

def clusters2rep(clusters):
    rep = {}
    for key in clusters:
        rep[key] = (major_vote(clusters[key]))
    return rep


def quality(constructed_rep, reference, type='sum'):
    if type == 'sum':
        ref = set(reference)
        con = set(constructed_rep)

        ref_con_intersection = ref.intersection(con)

        precision = 1.0 * len(ref_con_intersection) / len(con)
        sensitivity = 1.0 * len(ref_con_intersection) / len(ref)
        return (precision + sensitivity) / 2

    elif type == 'F1':
        return precision_sensetivity_F1(constructed_rep, reference)

    elif type == 'mult':
        ref = set(reference)
        con = set(constructed_rep)

        ref_con_intersection = ref.intersection(con)

        precision = 1.0 * len(ref_con_intersection) / len(con)
        sensitivity = 1.0 * len(ref_con_intersection) / len(ref)
        return precision * sensitivity


def clusters_classification(clusters, reference, constructed_rep):
    print "Repertoire construction started"
    #constructed_rep = clusters2rep(clusters)
    print "END"
    quality_0 = quality(constructed_rep.values(), reference)
    res = {}
    new_num = max([int(key) for key in clusters.keys()]) + 1
    for i, key in enumerate(clusters):
        if i%1000 == 0:
            print i
        if not second_vote(clusters[key]):
            res[key] = 0
            continue
        first_part, second_part = split_by_2nd_vote(clusters[key])
        temp_clusters = constructed_rep.copy()
        temp_clusters[key] = major_vote(first_part)
        temp_clusters[new_num] = major_vote(second_part)
        curr_quality = precision_sensetivity_F1(temp_clusters.values(), reference)
        if curr_quality > quality_0:
            quality_0 = curr_quality
            res[key] = 1
        elif curr_quality < quality_0:
            res[key] = -1
        else:
            res[key] = 0
    return res


def simple_clusters_classification(clusters, reference, constructed_rep):
    ref = set(reference)
    res = {}
    for i, key in enumerate(clusters):

        if i % 100 == 0:
            print i

        if not second_vote(clusters[key]):
            res[key] = 0
            continue

        first_part, second_part = split_by_2nd_vote(clusters[key])
        first_cons, second_cons = major_vote(first_part), major_vote(second_part)
        cluster_major = major_vote(clusters[key])

        if cluster_major in ref:
            if ((first_cons in ref) and (second_cons not in ref)) or ((second_cons in ref) and (first_cons not in ref)):
                res[key] = -1

            elif ((first_cons in ref) and (second_cons in ref)):
                res[key] = 1

            elif ((first_cons not in ref) and (second_cons not in ref)):
                res[key] = -1
        else:
            if ((first_cons in ref) and (second_cons not in ref)) or ((second_cons in ref) and (first_cons not in ref)):
                res[key] = 1

            elif ((first_cons in ref) and (second_cons in ref)):
                res[key] = 1

            elif ((first_cons not in ref) and (second_cons not in ref)):
                res[key] = -1
    return res


def clusters_filtering(clusters, threshold=5):
    filtered_clusters = {}
    for key in clusters:
        if len(clusters[key]) > threshold:
            filtered_clusters[key] = clusters[key]
    return filtered_clusters