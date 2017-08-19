import pandas as pd 
import numpy as np
import re

from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, train_test_split, RandomizedSearchCV, KFold
from sklearn.preprocessing import StandardScaler
import xgboost as xgb 
from sklearn.pipeline import Pipeline
import func_tools
from collections import Counter
import pickle

print 'Reading files with FLU and AGE datasets:' 

# full data for FLU[21, 22, 23, 27] datasets 
flu_data = pd.read_csv('../data/FLU.csv')
print 'FLU.csv is loaded'
# full data for AGE datasets 
age_data = pd.read_csv('../data/AGE.csv')
print 'AGE.csv is loaded\n'

print 'filter clusters with small size:'
print 'threshold: 5'
# filter clusters by size 
valid_flu = flu_data[flu_data['size'] >= 5]
valid_age = age_data[age_data['size'] >= 5]
print 'Filtered\n'

# and concat datasets 
compilation = pd.concat([valid_age, valid_flu], ignore_index=True)


# find IG and AGE datasets
reg = re.compile('.*_IG_.*')
ig_datasets = filter(reg.match, list(Counter(compilation['dataset'])))

print 'datasets for compilation:'
for dataset in ig_datasets:
    print dataset

ig_datasets = compilation[compilation['dataset'].isin(ig_datasets)]

reg = re.compile('AGE.*')
age_datasets = filter(reg.match, list(Counter(compilation['dataset'])))

print 'datasets for compilation:'
for dataset in age_datasets:
    print dataset

age_datasets = compilation[compilation['dataset'].isin(age_datasets)]


# create feature with absolute second value vote
ig_datasets['second_vote_abs1'] = ig_datasets['value1'] * ig_datasets['size']
age_datasets['second_vote_abs1'] = age_datasets['value1'] * age_datasets['size']

ig_datasets['second_vote_abs2'] = ig_datasets['value2'] * ig_datasets['size']
age_datasets['second_vote_abs2'] = age_datasets['value2'] * age_datasets['size']

ig_datasets['second_vote_abs3'] = ig_datasets['value3'] * ig_datasets['size']
age_datasets['second_vote_abs3'] = age_datasets['value3'] * age_datasets['size']


# construct mixed dataset
temp_df = pd.concat([ig_datasets[['second_vote_abs1', 'second_vote_abs2', 'second_vote_abs3', 'size']],
                    age_datasets[['second_vote_abs1', 'second_vote_abs2', 'second_vote_abs3', 'size']]])

scl = StandardScaler()
scl.fit(temp_df)


# use only two features
X_ig = pd.DataFrame(scl.transform(ig_datasets[['second_vote_abs1', 'second_vote_abs2', 'second_vote_abs3', 'size']]), 
                                  index=ig_datasets.index)

y_ig = pd.Series(ig_datasets['quality_imp'], index=ig_datasets.index).map(lambda x: 1 if x==1 else 0)

X_age = pd.DataFrame(scl.transform(age_datasets[['second_vote_abs1', 'second_vote_abs2', 'second_vote_abs3', 'size']]),
                                   index=age_datasets.index)

y_age  = pd.Series(age_datasets['quality_imp'], index=age_datasets.index).map(lambda x: 1 if x==1 else 0)


# sizes for samples weights
sizes_ig = ig_datasets['size']
sizes_age = age_datasets['size']

# mixed datasets
size_mixed = pd.concat([sizes_age, sizes_ig], ignore_index=True)
X_mixed = pd.concat([X_age, X_ig], ignore_index=True)
y_mixed = pd.concat([y_age, y_ig], ignore_index=True)

temp = pd.concat([X_mixed, y_mixed, size_mixed], axis=1).sample(frac=1)

X_mixed = temp[[0,1,2,3]]
y_mixed = temp['quality_imp'].map(lambda x: x if x == 1 else 0)

print 'Stacking model builder started'
models = func_tools.logreg_xgboost_stack(X_mixed, y_mixed, temp['size'])

first_lvl = func_tools.First_lvl_stacking(models)

pp = Pipeline(steps=[('scaler', scl), ('first_level', first_lvl), ('second_lvl', models['second_lvl'][0])])

pickle.dump(pp, open('../data/models/pipe_model', 'wb'))
