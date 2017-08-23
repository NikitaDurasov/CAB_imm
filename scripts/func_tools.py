import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter, defaultdict, OrderedDict
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, train_test_split, RandomizedSearchCV, KFold
import re
import itertools
import operator
import xgboost as xgb

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
    return map(lambda x: str(x.seq), SeqIO.parse(filename, "fasta"))

def read_rcm(filename):
    """
    Read .rcm file into dict obj "read_id" -> "cluster_id"
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
    :return: dict odj with (position in cluster) -> (Counter obj) or 0 if there is no any variations
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
    For every position in cluster function calculates second vote value - frequency of second the most frequent letter
    in position
    :param cluster: list of reads
    :return: dict odj with (position in cluster) -> (second vote value)
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


def second_vote_letter(cluster, position):
    if not np.isnan(position):
        variance = cluster_variety(cluster)[position]
        if variance:
            temp = sorted(variance.items(), key=lambda x: x[1])
            return temp[0][0]
        else:
            return ' '
    else:
        return ' '


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

def to_fa(filename, clusters):
    rep = clusters2rep(clusters)
    sizes = clusters_size_dict(clusters)
    res = {}

    for key in clusters:
        res[key] = [rep[key], sizes[key]]

    f = open(filename, 'w')

    for key in res:
        f.write('>cluster___' + key + '___size___' + str( res[key][1]) + '\n')
        temp_str = res[key][0]

        while temp_str:
            f.write(temp_str[:60] + '\n')
            temp_str = temp_str[60:]

def precision_sensitivity_F1(constructed, reference):
    """
    Calculates F1 for precision and sensitivity for constructed and reference repertoire
    :param constructed: constructed repertoire - list of strings
    :param reference: reference repertoire - list of strings
    :return: double value F1
    """
    ref = set(reference)
    con = set(constructed)

    ref_con_intersection = ref.intersection(con)

    precision = 1.0 * len(ref_con_intersection) / len(con)
    sensitivity = 1.0 * len(ref_con_intersection) / len(ref)

    return 2.0 * precision * sensitivity / (precision + sensitivity)


def split_cluster(cluster, position, threshold=5):
    """
    Splits cluster into two clusters by position using second vote value
    :param cluster: list of reads
    :param position: int value
    :return: list of parts of divided clusters; if new cluster cardinality is lower than 5, then this splitted part
    turned to empty list
    """
    position_variation = cluster_variety(cluster)[position]
    if not position_variation and len(cluster) >= 5:
        return (cluster, [])
    elif not position_variation and not len(cluster) >=5:
        return ([],[])
    letter_for_split = find_max_dict(position_variation)
    first_splitted_part = []
    second_splitted_part = []

    for item in cluster:
        if item[position] == letter_for_split:
            first_splitted_part.append(item)
        else:
            second_splitted_part.append(item)

    return check_size(first_splitted_part, second_splitted_part, threshold=threshold)


def find_max_dict(dictionary):
    """
    Find maximal element in dictionary in values
    :param dictionary: dict with numeric value
    """
    return max(dictionary.iteritems(), key=operator.itemgetter(1))[0] if dictionary else -1


def check_size(first_splitted_part, second_splitted_part, threshold=5):
    if len(first_splitted_part) < threshold and len(second_splitted_part) < threshold:
        return ([], [])

    elif len(first_splitted_part) < threshold and len(second_splitted_part) >= threshold:
        return (second_splitted_part, [])

    elif len(first_splitted_part) >= threshold and len(second_splitted_part) < threshold:
        return (first_splitted_part, [])

    else:
        return (first_splitted_part, second_splitted_part)

def split_by_2nd_vote(cluster, threshold=5):
    max_2nd_vote_pos = find_max_dict(second_vote(cluster))
    return split_cluster(cluster, max_2nd_vote_pos, threshold=threshold) if max_2nd_vote_pos != -1 else check_size(cluster, [], threshold=threshold)

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
        return precision_sensitivity_F1(constructed_rep, reference)

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
        curr_quality = precision_sensitivity_F1(temp_clusters.values(), reference)
        if curr_quality > quality_0:
            quality_0 = curr_quality
            res[key] = 1
        elif curr_quality < quality_0:
            res[key] = -1
        else:
            res[key] = 0
    return res


def simple_clusters_classification(clusters, reference):
    ref = set(reference)
    res = {}
    for i, key in enumerate(clusters):

        if i % 1000 == 0:
            print i, ' reads checked'

        if not second_vote(clusters[key]):
            res[key] = [-1, '?', '?', '?']
            continue

        first_part, second_part = split_by_2nd_vote(clusters[key], threshold=5)
        cluster_major = major_vote(clusters[key])

        #TO DO: REWRITE
        if len(first_part)*len(second_part):

            first_cons, second_cons = major_vote(first_part), major_vote(second_part)

            if cluster_major in ref:
                if ((first_cons in ref) and (second_cons not in ref)) or ((second_cons in ref) and (first_cons not in ref)):
                    res[key] = [-1, '+', '+', '-']

                elif first_cons in ref and second_cons in ref:
                    res[key] = [1, '+', '+', '+']

                elif first_cons not in ref and second_cons not in ref:
                    res[key] = [-1, '+', '-', '-']
            else:
                if ((first_cons in ref) and (second_cons not in ref)) or ((second_cons in ref) and (first_cons not in ref)):
                    res[key] = [-1, '-', '+', '-']

                elif first_cons in ref and second_cons in ref:
                    res[key] = [1, '-', '+', '+']

                elif first_cons not in ref and second_cons not in ref:
                    res[key] = [-1, '-', '-', '-']
        else:
            if len(first_part):
                if cluster_major in reference:
                    if major_vote(first_part) in reference:
                        res[key] = [-1, '+', '+', '?']
                    else:
                        res[key] = [-1, '+', '-', '?']
                else:
                    if major_vote(first_part) in reference:
                        res[key] = [-1, '-', '+', '?']
                    else:
                        res[key] = [-1, '-', '-', '?']

            elif len(second_part):
                if cluster_major in reference:
                    if major_vote(second_part) in reference:
                        res[key] = [-1, '+', '?', '+']
                    else:
                        res[key] = [-1, '+', '?', '-']
                else:
                    if major_vote(second_part) in reference:
                        res[key] = [-1, '-', '?', '+']
                    else:
                        res[key] = [-1, '-', '?', '-']
            else:
                if cluster_major in reference:
                    res[key] = [-1, '+', '?', '?']
                else:
                    res[key] = [-1, '-', '?', '?']
        #END


    return res

def clusters_filtering(clusters, threshold=5):
    filtered_clusters = {}

    for key in clusters:
        if len(clusters[key]) > threshold:
            filtered_clusters[key] = clusters[key]

    return filtered_clusters


def n_second_vote(second_vote, n=0):
    res  = {}
    for key in second_vote:
        if len(second_vote[key].values()) > n:
            temp = sorted(second_vote[key].items(), key=lambda x: x[1], reverse=True)
            res[key] = temp[n]
        else:
            res[key] = (np.nan, 0)
    return res


def find_context(repertoire, max_second_vote, n=2):
    res = {}
    for key in repertoire:
        position = max_second_vote[key][0]
        if not np.isnan(position):
            sequence = " "*n + repertoire[key] + " "*n
            res[key] = sequence[position:position+2*n+1]
        else:
            res[key] = " "*(2*n+1)
    return res


def build_df(input_reads, rcm_file, rcm_reference=None, fa_reference=None,
        classification='simple', threshold=10):
    print 'build_df started'

    id_dict = id_to_read(input_reads)

    if fa_reference:
        rep = read_repertoire(fa_reference)

    if rcm_reference:
        reference_rcm = read_rcm(rcm_reference)

    igrec_rcm = read_rcm(rcm_file)

    igrec_clusters = construct_clusters(igrec_rcm, id_dict)
    print 'All clusters: ', len(igrec_clusters)
    igrec_clusters = clusters_filtering(igrec_clusters, threshold=threshold)

    igrec_rep = clusters2rep(igrec_clusters)
    igrec_res = second_votes(clusters_filtering(igrec_clusters, threshold=0))

    print 'calculation step'

    max_final_second_vote = n_second_vote(igrec_res)
    max_2nd_final_second_vote = n_second_vote(igrec_res, n=1)
    max_3nd_final_second_vote = n_second_vote(igrec_res, n=2)
    sizes = clusters_size_dict(igrec_clusters)
    #context = find_context(clusters_filtering(igrec_rep), max_final_second_vote)

    #second_vote_std = {}
    #for key in igrec_res:
    #    second_vote_std[int(key)] = np.std(igrec_res[key].values())

    print 'parsing step'

    #pos1 = {int(k): v[0] for k, v in max_final_second_vote.items()}
    value1 = {int(k): v[1] for k, v in max_final_second_vote.items()}
    #pos2 = {int(k): v[0] for k, v in max_2nd_final_second_vote.items()}
    value2 = {int(k): v[1] for k, v in max_2nd_final_second_vote.items()}
    #pos3 = {int(k): v[0] for k, v in max_3nd_final_second_vote.items()}
    value3 = {int(k): v[1] for k, v in max_3nd_final_second_vote.items()}
    sizes = {int(k): v for k, v in sizes.items()}

    #print 'context step'

    #context1 = {}
    #context2 = {}
    #context3 = {}
    #context4 = {}
    #context5 = {}

    #for key in context:
    #    temp_arr = list(context[key])
    #    context1[int(key)] = temp_arr[0]
    #    context2[int(key)] = temp_arr[1]
    #    context3[int(key)] = temp_arr[2]
    #    context4[int(key)] = temp_arr[3]
    #    context5[int(key)] = temp_arr[4]

    #print 'mutation step'

    #mutated_letter = {}
    #for key in igrec_res:
    #    mutated_letter[int(key)] = second_vote_letter(igrec_clusters[key], max_final_second_vote[key][0])

    df = pd.DataFrame({#'pos1': pos1, 
                       'value1': value1,
                       #'pos2': pos2, 
                       'value2': value2,
                       #'pos3': pos3, 
                       'value3': value3,
                       #'context1': context1,
                       #'context2': context2,
                       #'context3': context3,
                       #'context4': context4,
                       #'context5': context5,
                       #'mutated_letter': mutated_letter,
                       'size': sizes
                       #'second_vote_std': second_vote_std})
                        })
    df.index = df.index.map(int)
    final_df = df

    if fa_reference:

        print 'answers DataFrame'

        if classification == 'simple':
            ans_dict = simple_clusters_classification(igrec_clusters, rep)
            ans_df = build_ans_df(ans_dict)

            final_df = pd.concat([final_df, ans_df], axis=1)

        elif classification == 'colormap':
            ans_dict, _ = reference_classification(igrec_clusters, reference_rcm, rep, id_dict)
            ans_df = pd.DataFrame({'quality_imp': ans_dict})
            ans_df.index = ans_df.index.map(int)

            final_df = pd.concat([final_df, ans_df], axis=1)

    print 'building succeeded'
    return final_df


def build_df_preload(igrec_clusters, igrec_rep, igrec_res):

    print 'build_df started'
    print 'calculation step'

    max_final_second_vote = n_second_vote(igrec_res)
    max_2nd_final_second_vote = n_second_vote(igrec_res, n=1)
    max_3nd_final_second_vote = n_second_vote(igrec_res, n=2)
    sizes = clusters_size_dict(clusters_filtering(igrec_clusters))
    context = find_context(clusters_filtering(igrec_rep), max_final_second_vote)

    second_vote_std = {}
    for key in igrec_res:
        second_vote_std[key] = np.std(igrec_res[key].values())

    print 'parsing step'

    pos1 = {k: v[0] for k,v in max_final_second_vote.items()}
    value1 = {k: v[1] for k,v in max_final_second_vote.items()}
    pos2 = {k: v[0] for k,v in max_2nd_final_second_vote.items()}
    value2 = {k: v[1] for k,v in max_2nd_final_second_vote.items()}
    pos3 = {k: v[0] for k,v in max_3nd_final_second_vote.items()}
    value3 = {k: v[1] for k,v in max_3nd_final_second_vote.items()}

    print 'context step'

    context1 = {}
    context2 = {}
    context3 = {}
    context4 = {}
    context5 = {}

    for key in context:
        temp_arr = list(context[key])
        context1[key] = temp_arr[0]
        context2[key] = temp_arr[1]
        context3[key] = temp_arr[2]
        context4[key] = temp_arr[3]
        context5[key] = temp_arr[4]

    print 'mutation step'

    mutated_letter = {}
    for key in igrec_res:
        mutated_letter[key] = second_vote_letter(igrec_clusters[key], max_final_second_vote[key][0])

    df = pd.DataFrame({'pos1':pos1, 'value1':value1,
                        'pos2': pos2, 'value2': value2,
                        'pos3': pos3, 'value3': value3,
                        'context1': context1,
                        'context2': context2,
                        'context3': context3,
                        'context4': context4,
                        'context5': context5,
                        'mutated_letter': mutated_letter,
                        'size': sizes,
                        'second_vote_std':second_vote_std})

    df.index = df.index.map(int)

    print 'building succeeded'

    return df

def build_ans_df(answer_dict):

    class_label = {}
    origin_cluster = {}
    first_cluster = {}
    second_cluster = {}

    for key in answer_dict:
        temp_list = answer_dict[key]
        class_label[key] = temp_list[0]
        origin_cluster[key] = temp_list[1]
        first_cluster[key] = temp_list[2]
        second_cluster[key] = temp_list[3]

    ans_df = pd.DataFrame({'quality_imp':class_label,
                           'origin_cluster': origin_cluster,
                           'first_cluster': first_cluster,
                           'second_cluster': second_cluster})

    ans_df.index = ans_df.index.map(int)

    return ans_df

def clusters_splitting(clusters, target_clusters, threshold=5):

    temp_clusters = clusters.copy()
    clusters_id = map(lambda x: int(x), clusters.keys())
    max_id = max(clusters_id) + 1

    print "LENGHT OF TARGET_CLUSTERS: ", len(target_clusters)

    for cluster in target_clusters:
        print 'Splitting cluster: ', cluster, ' with length ', len(temp_clusters[str(cluster)])
        first_part, second_part = split_by_2nd_vote(clusters[str(cluster)], threshold=0)

        if len(first_part) > threshold and len(second_part) > threshold:

            del temp_clusters[str(cluster)]

            if first_part and not second_part:
                temp_clusters[str(max_id)] = first_part
                print 'New cluster ', max_id, ' with length ', len(first_part)
                print '\n'
                max_id += 1

            elif second_part and not first_part:
                temp_clusters[str(max_id)] = second_part
                print 'New cluster ', max_id, ' with length ', len(second_part)
                print '\n'
                max_id += 1

            elif second_part and second_part:
                temp_clusters[str(max_id)] = first_part
                temp_clusters[str(max_id + 1)] = second_part
                print 'New clusters ', max_id, ' and ', max_id + 1,  ' with length ', len(first_part), ' and ', len(second_part)
                max_id += 2
                print '\n'
            else:
                print 'Both parts are empty\n'

        else:
            print "One/Both of splitted clusters is/are too small\n"

    return temp_clusters

def to_rcm(filename, clusters, id_dict):

    #read2id = {v:k for k, v in id_dict.items()}
    #output = open(filename, 'w')

    #for cluster in clusters:
    #    for read in clusters[cluster]:
    #        output.write(read2id[read] + '\t' + cluster + '\n')

    output = open(filename, 'w')

    read2id = defaultdict(lambda : [])
    for read_id in id_dict:
        read2id[id_dict[read_id]].append(read_id)

    for cluster in clusters:
        for read in clusters[cluster]:
            output.write(read2id[read][0] + '\t' + cluster + '\n')
            del read2id[read][0]

    output.close()


class First_lvl_stacking:

    def __init__(self, models):
        self.models = models

    def fit(self, x):
        pass

    def transform(self, x):
        logreg_pred = self.models['first_lvl'][0].predict(x)
        xgb_pred = self.models['first_lvl'][1].predict(x)
        return np.matrix([logreg_pred, xgb_pred]).T


def logreg_xgboost_stack(X, y, sample_weight=None):

        logreg_cls = LogisticRegression()
        xgb_cls = xgb.XGBClassifier()

        X_train_first, X_train_second, y_train_first, y_train_second = train_test_split(X, y, stratify=y,
        test_size=0.2, train_size=0.8)

        kf = KFold(n_splits=5, shuffle=True)

        ans_logreg = pd.Series()
        ans_xgb = pd.Series()

        best_param1 = [{}, 0]
        best_param2 = [{}, 0]

        gs1 = GridSearchCV(estimator=logreg_cls, param_grid={'C': np.linspace(1e-4, 1, 10)},
                               scoring='roc_auc', fit_params={'sample_weight': sample_weight.iloc[X_train_first.index]}, n_jobs=1)
        gs2 = RandomizedSearchCV(estimator = xgb_cls,
                                 scoring = 'roc_auc',
                                 param_distributions = {'subsample': np.linspace(0.5, 0.9, 3),
                                                        'max_depth': range(1,6,2), 
                                                        'n_estimators': range(20, 70, 10), 
                                                        'min_child_weight': range(2,8,3), 
                                                        'learning_rate': np.linspace(1e-4, 0.1, 5)}, 
                                 cv=3, n_jobs=1,
                                 fit_params={'sample_weight': sample_weight.iloc[X_train_first.index]}) 

        gs1.fit(X_train_first, y_train_first)

        gs2.fit(X_train_first, y_train_first)

        pred_feat_logreg = pd.Series(gs1.best_estimator_.predict(X_train_second), index=X_train_second.index)
        pred_feat_xgb = pd.Series(gs2.best_estimator_.predict(X_train_second), index=X_train_second.index)

        logreg_cls_second = LogisticRegression()
        logreg_cls_second.fit(pd.concat([pred_feat_logreg, pred_feat_xgb], axis=1), y_train_second)

        return {'first_lvl': [logreg_cls_first, xgb_cls_first], 'second_lvl': [logreg_cls_second]}


def reference_classification(clusters, reference_rcm, reference, id_dict):

    res = {}
    read_to_id = {v:k for k,v in id_dict.items()}
    clusters_perc = {}

    for cluster in clusters:
        cl_size = len(clusters[cluster])

        if cl_size < 10:
            res[cluster] = 0

        else:
            cluster_parts = []
            for read in clusters[cluster]:
                read_id = read_to_id[read]
                ref_cluster_id = reference_rcm[read_id]
                cluster_parts.append(ref_cluster_id)

            freq = Counter(cluster_parts)
            freq = sorted(freq.items(), key=lambda x: x[1], reverse=True)
            relative_freq = [(k, 1.0 * v / cl_size) for k,v in freq]

            if len(freq) == 1:
                res[cluster] = 0

            elif (freq[0][1] >= 5 and freq[1][1] >= 5):

                clusters_perc[cluster] = [relative_freq[0][1],
                                           relative_freq[1][1]]

                first_part = []
                second_part = []

                for read in clusters[cluster]:
                    read_id = read_to_id[read]

                    if reference_rcm[read_id] == freq[0][0]:
                        first_part.append(read)

                    if reference_rcm[read_id] == freq[1][0]:
                        second_part.append(read)

                if (major_vote(first_part) in reference and
                    major_vote(second_part) in reference):

                    res[cluster] = 1

                else:
                    res[cluster] = 0

            else:
                res[cluster] = 0
    return (res, clusters_perc)

def find_unrecognized_clusters(reference_repertoire, constructed_repertoire):

    res = []
    rep, cons_rep = set(reference_repertoire.values()), set(constructed_repertoire.values())
    clusters_intersection = set.intersection(rep, cons_rep)

    for cluster in reference_repertoire:
        if reference_repertoire[cluster] not in clusters_intersection:
            res.append(cluster)

    return res

def compare_sizes(first_rcm_file, second_rcm_file, input_reads):

    id_dict = id_to_read(input_reads)

    first_rcm = read_rcm(first_rcm_file)
    second_rcm = read_rcm(second_rcm_file)

    first_clusters = construct_clusters(first_rcm, id_dict)
    second_clusters = construct_clusters(second_rcm, id_dict)

    first_rep = clusters2rep(first_clusters)
    second_rep = clusters2rep(second_clusters)

    first_reads = set(first_rep.values())
    second_reads = set(second_rep.values())

    reads_intersection = first_reads.intersection(second_reads)

    #FINISH

def sens_prec_plot(igrec_json, test_json, axe):
        axe.plot(test_json['reference_based']['__data_precision'], test_json['reference_based']['__data_sensitivity'], label='test', c='r')
        axe.plot(igrec_json['reference_based']['__data_precision'], igrec_json['reference_based']['__data_sensitivity'], label='igrec', c='b')

        axe.scatter(test_json['reference_based']['__data_precision'], test_json['reference_based']['__data_sensitivity'], s=5, c='r')
        axe.scatter(igrec_json['reference_based']['__data_precision'], igrec_json['reference_based']['__data_sensitivity'], s=5, c='b')

        axe.scatter(test_json['reference_based']['__data_precision'][4], test_json['reference_based']['__data_sensitivity'][4], s=20, c='k')
        axe.scatter(igrec_json['reference_based']['__data_precision'][4], igrec_json['reference_based']['__data_sensitivity'][4], s=20, c='k')


def res_cons(input_reads, ref_rcm_file, cons_rcm_file):
        id_dict = func_tools.id_to_read(input_reads)

        ref_rcm = func_tools.read_rcm(ref_rcm_file)
        cons_rcm = func_tools.read_rcm(cons_rcm_file)

        ref_clusters = func_tools.construct_clusters(ref_rcm, id_dict)
        cons_clusters = func_tools.construct_clusters(cons_rcm, id_dict)

        rep = func_tools.clusters2rep(ref_clusters)
        cons_rep = func_tools.clusters2rep(cons_clusters)

        res = func_tools.find_unrecognized_clusters(rep, cons_rep)

        return (res, ref_clusters, cons_clusters)

def unrecognized_clusters_sizes(input_reads, ref_rcm_file, cons_rcm_file):
    res = res_cons(input_reads, ref_rcm_file, cons_rcm_file)

    r = set(res[0])
    ref_clusters = res[1]
    cons_clusters = res[2]

    missed_clusters = {cluster:ref_clusters[cluster] for cluster in r}
    missed_sizes = func_tools.clusters_size_dict(missed_clusters)
    all_sizes = func_tools.clusters_size_dict(ref_clusters)
    missed_freq = defaultdict(lambda : 0, Counter(missed_sizes.values()))
    all_freq = defaultdict(lambda : 0, Counter(all_sizes.values()))

    perc_sizes_hist = []
    for i in range(5, 200):
        if all_freq[i]:
            sizes_hist.append(1.0 * missed_freq[i] / all_freq[i])
        else:
            sizes_hist.append(0)

    sizes_hist = []
    for i in range(5, 200):
        if all_freq[i]:
            sizes_hist.append(missed_freq[i])
        else:
            sizes_hist.append(0)

    return_dict = {}
    return_dict['perc_sizes_hist'] = perc_sizes_hist
    return_dict['sizes_hist'] = sizes_hist
    return_dict['reference_clusters'] = ref_clusters
    return_dict['constructed_clusters'] = cons_clusters 
    return_dict['unrecognized_clusters'] = r

    return return_dict
