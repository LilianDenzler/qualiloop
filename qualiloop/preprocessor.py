#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Preprocessor
   File:       preprocessor.py
   
   Version:    V1.1
   Date:       17.03.21
   Function:   
   
   Copyright:  (c) UCL, Prof. Andrew C. R. Martin 1994-2021
   Author:     Lilian Denzler
   Address:    Biomolecular Structure & Modelling Unit,
               Department of Biochemistry & Molecular Biology,
               University College,
               Gower Street,
               London.
               WC1E 6BT.
   EMail:      andrew@bioinf.org.uk, lilian.denzler.17@ucl.ac.uk
               
**************************************************************************
   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC
   The code may be modified as required, but any modifications must be
   documented so that the person responsible can be identified. If someone
   else breaks this code, I don't want to be blamed for code that does not
   work! 
   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.
**************************************************************************
   Description:
   ============
   The program contains all functions needed for different forms of data
   preprocessing and feature encoding.
   All of these differnet pipelines will be used to assess which types of 
   preprocessing works best with which model. 


**************************************************************************
   Usage:
   ======
   This library is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   21.01.21 Original
*************************************************************************/"""
import sys
import os
import pandas as pd
import numpy as np
import joblib
import time
import warnings
warnings.filterwarnings("ignore", message="../src/learner.cc:1061: Starting in XGBoost 1.3.0, the default evaluation metric used with the objective 'multi:softprob' was changed from 'merror' to 'mlogloss'. Explicitly set eval_metric if you'd like to restore the old behavior.")

import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import RepeatedStratifiedKFold

from statsmodels.stats.outliers_influence import variance_inflation_factor
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import make_scorer
import xgboost as xgb
import collections

from pathlib import Path
import numpy as np
from sklearn.decomposition import TruncatedSVD
from sklearn.decomposition import PCA, IncrementalPCA
import collections
from sklearn import manifold, datasets
from functools import partial
from sklearn.feature_selection import VarianceThreshold
from sklearn.feature_selection import SelectKBest
#from sklearn.feature_selection import RFECV
from yellowbrick.model_selection import RFECV 

from sklearn.svm import SVR
from sklearn import preprocessing
from sklearn.feature_selection import SelectFromModel
import featuretools as ft

from qualiloop import physiochem_4d_aa_encoder

import glob

#------------------------------
#LOCATING MATRIX FILE
def get_file_paths(matrix):
  matrices_path=glob.glob('**/data/**/matrices_blosum_pam/'+matrix, recursive=True)
  #print('**/matrices_blosum_pam/'+matrix)
  if len(matrices_path)>1:
    print(matrices_path)
    raise Warning('more than one matrices_blosum file was found')
  elif len(matrices_path)==0:
    print("matrix file is not in data!!!!")
    matrices_path=glob.glob('**/matrices_blosum_pam/'+matrix, recursive=True)
    raise Warning("matrix file is not in data!!!!")
    if len(rangecontacts_path)==0:
      raise ValueError('no matrix file was found')
  return matrices_path[0]

#-------------------------------------------------------------------------------------------------------------
#DATA SCALING AND CLEANING
def scaling(df, scaler_name, run_dir):

  if scaler_name == 'standard':
    scaler = StandardScaler()
  elif scaler_name =='minmax':
    scaler = MinMaxScaler()
  elif scaler_name =='robust':
    scaler= RobustScaler(with_centering=True,with_scaling=True,quantile_range=(25.0, 75.0),copy=True)

  columns_to_scale=list(df.select_dtypes(include=[np.number]).columns.values)
  columns_to_scale_feat=[col for col in columns_to_scale if "local_CA" not in col and "local_AA" not in col and "global_CA" not in col and "global_AA"  not in col]
  columns_to_scale=columns_to_scale_feat
  #select all numeric columns, that have more than two values ( i.e they are not binary)
  for i in columns_to_scale:
    if len(np.unique(df[i])) <= 2:
      columns_to_scale.remove(i)

  df[columns_to_scale] = scaler.fit_transform(df[columns_to_scale].to_numpy())
  print("SCALED")
  with open(os.path.join(run_dir,"info.txt"), 'a') as file:
      file.write("Scaling numerical features using scaler: {}\n".format(scaler_name))

  return df


def remove_outliers(df, run_dir):
  columns_to_scale=list(df.select_dtypes(include=[np.number]).columns.values)
  columns_to_scale_new=[col for col in columns_to_scale if "local_CA" not in col and "local_AA" not in col and "global_CA" not in col and "global_AA"  not in col]
  outlier_indexes_list=[]
  df.reset_index()
  for feature_name in columns_to_scale_new: 
    feature_column = df[feature_name]
    Q1 = np.percentile(feature_column, 25.) ###changed from 25/75 as this gives 619 outliers
    Q3 = np.percentile(feature_column, 75.)

    #An outlier is defined as being any point of data that lies over 1.5 
    #IQRs(interquartile range) below the first quartile (Q1) or above the third quartile (Q3)in a data set.
    #High = (Q3) + 1.5 IQR
    #Low = (Q1) – 1.5 IQR

    # Use the interquartile range to calculate an outlier step (1.5 times the interquartile range)
    step = (Q3-Q1)*1.5
    outlier_indexs = feature_column[~(feature_column >= (Q1 - step)) & (feature_column <= (Q3 + step))].index.tolist()
    print(feature_name, len(outlier_indexs))
    outlier_indexes_list.extend(outlier_indexs)
    
  outlier_indexes_list=np.array(list(set(tuple(outlier_indexes_list))))
  with open(os.path.join(run_dir,"info.txt"), 'a') as file:
      #file.write("Outlier dataframe:\n {}".format(df.iloc[outlier_indexes_list.tolist()]))
      file.write("Number of utliers: {}\n".format(len(outlier_indexes_list)))
  new_dataset_df=df.drop(outlier_indexes_list).reset_index(drop = True)
 

  return new_dataset_df



#-----------------------------------------------------------------------------------------------------------------------------
#AMINO ACID ENCODING
def aa_encoding_all(full_dataset_df, non_num_cols, run_dir):
  #['seq', 'tip_res', 'tip', 'res', 'atom_name', 'chain', 'pos', 'critical_seq']
  for column in non_num_cols:
    print(column)
    col_arr=full_dataset_df[column].to_numpy()
    one_hot_col=pd.DataFrame()
    blosum_62_col=pd.DataFrame()
    NLF_col=pd.DataFrame()
    physiochem_4d_col=pd.DataFrame()  
    for i in col_arr:
      if i ==np.nan or i==None:
        one_hot_col=pd.concat([one_hot_col, [np.nan]*len(one_hot_col.columns)], axis=0) #if only one dimensional (i.e. only one residue is encoded) no dimensionality reduction
        blosum_62_col=pd.concat([blosum_62_col, [np.nan]*len(blosum_62_col)],axis=0)
        NLF_col=pd.concat([NLF_col, [np.nan]*len(NLF_col)],axis=0)
        physiochem_4d_col=pd.concat([physiochem_4d_col, [np.nan]*len(physiochem_4d_col)],axis=0)
        continue
      one_hot_df=one_hot(i)
      blosum_62_df=blosum62(i)
      NLF_df=NLF_encode(i)
      physiochem_4d_df=physiochem_4d(i)
      if physiochem_4d_df.empty or NLF_df.empty or blosum_62_df.empty or one_hot_df.empty:
        input("EMPTY_ ROW")
      one_hot_col=pd.concat([one_hot_col, one_hot_df], axis=0) #if only one dimensional (i.e. only one residue is encoded) no dimensionality reduction
     
      blosum_62_col.columns.duplicated()
      blosum_62_col.T.columns.duplicated()
      blosum_62_col=pd.concat([blosum_62_col, blosum_62_df.T],axis=0)
      NLF_col=pd.concat([NLF_col, NLF_df],axis=0)
      physiochem_4d_col=pd.concat([physiochem_4d_col, physiochem_4d_df],axis=0)

    #one_hot_col=reduction_set(one_hot_df, column+"_one_hot") #Transpose, as we want to reduce the nr of columns (i.e. below 20)
    blosum_62_col=reduction_set(blosum_62_col, column+"_blosum_62", run_dir)
    NLF_col=reduction_set(NLF_col, column+"_NLF", run_dir)

    one_hot_col.columns=[column+"_one_hot_" for column in one_hot_col.columns]
    blosum_62_col.columns=[column+"_blosum_62_" for column in blosum_62_col.columns]
    NLF_col.columns=[column+"_NLF_col_" for column in NLF_col.columns]
    physiochem_4d_col.columns=[column+"_physiochem_4d_" for column in physiochem_4d_col.columns]
    
    print(one_hot_col)
    print(blosum_62_col)
    print(NLF_col)
    print(physiochem_4d_col)
    del full_dataset_df[column]
    print("full_df shape:",full_dataset_df.shape, "  one_hot shape:", one_hot_col.shape, " blosum_62 shape", blosum_62_col.shape, "NLF", NLF_col.shape, "physiochem_4d_df", physiochem_4d_col.shape)
    full_dataset_df.reset_index(drop=True, inplace=True)
    one_hot_col.reset_index(drop=True, inplace=True)
    blosum_62_col.reset_index(drop=True, inplace=True)
    NLF_col.reset_index(drop=True, inplace=True)
    physiochem_4d_col.reset_index(drop=True, inplace=True)
    one_hot_col=one_hot_col.iloc[:, 1:]
    full_dataset_df=pd.concat([full_dataset_df, one_hot_col], axis=1)
    full_dataset_df=pd.concat([full_dataset_df, blosum_62_col], axis=1)
    full_dataset_df=pd.concat([full_dataset_df, NLF_col], axis=1)
    full_dataset_df=pd.concat([full_dataset_df, physiochem_4d_col], axis=1)
    with open(os.path.join(run_dir,"info.txt"), 'a') as file:
      file.write("""Amino acid encoding of the feature: {} using One-Hot encoding, BLOSUM62 encoding, NLF encoding and physicochemical 
        4-d vector encoding.\n""".format(column))
  return full_dataset_df

def one_hot(seq):
    amino_acids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L','M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    o = list(set(amino_acids) - set(seq))
    s = pd.DataFrame(list(seq))    
    x = pd.DataFrame(np.zeros((len(seq),len(o)),dtype=int),columns=o)    
    a = s[0].str.get_dummies(sep=',')
    a = a.join(x)
    a = a.sort_index(axis=1)
    one_hot=a
    #e = a.values.flatten() #dont do if you want as df
    #one_hot=e
    if one_hot.empty ==True:
      input("EMPTY ONE_HOT_DF")
    return one_hot

def blosum_pam_score(res1, res2, matrix):
  filepath="data/matrices_blosum_pam/BLOSUM62"
  matrix=[]
  try:
    with open(filepath) as file:
      for line in file:
        if not line.startswith("#"):
          #print(line)
          line=line.split()
          matrix.append(line)
  except:
    filepath=get_file_paths(matrix)
    with open(filepath) as file:
      for line in file:
        if not line.startswith("#"):
          #print(line)
          line=line.split()
          matrix.append(line)
  columns=matrix[0]
  matrix=[line[1:] for line in matrix] 
  matrix=matrix[1:]
  matrix_df=pd.DataFrame(data=matrix, index=columns, columns=columns)
  score=matrix_df[res1][res2]
  return(score)


def blosum62(seq):
  filepath="data/matrices_blosum_pam/BLOSUM62"
  matrix=[]
  try:
    with open(filepath) as file:
      for line in file:
        if not line.startswith("#"):
          line=line.split()
          matrix.append(line)
  except:
    filepath=get_file_paths("BLOSUM62")
    with open(filepath) as file:
      for line in file:
        if not line.startswith("#"):
          line=line.split()
          matrix.append(line)

  columns=matrix[0]
  matrix=[line[1:] for line in matrix] 
  matrix=matrix[1:]
  matrix_df=pd.DataFrame(data=matrix, index=columns, columns=columns)
  blosum_62_df=pd.DataFrame()
  for i in seq:
    res_column=matrix_df[i]
    blosum_62_df=pd.concat([blosum_62_df,res_column], axis=0)
  if blosum_62_df.empty ==True:
    input("EMPTY BLOSUM_DF")
  return blosum_62_df


def NLF_encode(seq):
  #his method of encoding is detailed by Nanni and Lumini in their paper. It takes many physicochemical properties and transforms them using a Fisher Transform 
  #(similar to a PCA) creating a smaller set of features that can describe the amino acid just as well. There are 19 transformed features. 
  #nlf = pd.read_csv('https://raw.githubusercontent.com/dmnfarrell/epitopepredict/master/epitopepredict/mhcdata/NLF.csv',index_col=0)
  nlf = pd.read_csv('data/NLF.csv',index_col=0)
  x = pd.DataFrame([nlf[i] for i in seq]).reset_index(drop=True)  
  #e = x.values.flatten()
  if x.empty ==True:
    input("EMPTY NLF_DF")
  return x

def zscale(feature_directory,input_seq_directory):
    pass
def VHSE_scale(feature_directory,input_seq_directory):
    pass
def BLOSUM_50():
    pass

def physiochem_4d(seq):
  physiochem_4d_df=physiochem_4d_aa_encoder.run(seq)
  if physiochem_4d_df.empty ==True:
    input("EMPTY physiochem_DF")
  return physiochem_4d_df
  #Martin&Abhinandan 2010. the amino acids were encoded in 4d vectors for neural network input preparation
  #1. total number of side-chain atoms
  #2. number of side-chain atoms in shortest path from Calpha to most distal atom
  #3. eisenberg consensus hydrophobicity 
  #4. charge (histidine was assigned +0.5)


#-----------------------------------------------------------------------------------------------------------------------------------
#ENCODING TIP RESIDUE POSITION

def encode_tip(full_dataset_df,column):
  del full_dataset_df[column]
  return full_dataset_df
#-----------------------------------------------------------------------------------------------------------------------------------
#DIMENSIONALITY REDUCTION

def trunc_SVD(df):
  svd = TruncatedSVD(n_components=2)
  svd.fit(df)
  result = svd.transform(df)
  return result

def reduction_set(df,feature_name, run_dir):
  n_neighbors = 5
  n_components=3
  LLE = partial(manifold.LocallyLinearEmbedding,n_neighbors, n_components, eigen_solver='auto')         #https://scikit-learn.org/stable/modules/manifold.html#manifold
  methods = collections.OrderedDict()
  methods['LLE'] = LLE(method='standard')
  #methods['LTSA'] = LLE(method='ltsa')
  methods['Mod_LLE'] = LLE(method='modified', eigen_solver='dense')
  methods['Isomap'] = manifold.Isomap(n_neighbors, n_components)
  methods['MDS'] = manifold.MDS(n_components, max_iter=100, n_init=1)
  methods['SE'] = manifold.SpectralEmbedding(n_components=n_components,n_neighbors=n_neighbors)        #https://scikit-learn.org/stable/auto_examples/manifold/plot_compare_methods.html
  methods['tSNE'] = manifold.TSNE(n_components=n_components, init='pca',random_state=0)               #https://dmnfarrell.github.io/bioinformatics/mhclearning
  methods["PCA"]= IncrementalPCA(n_components=n_components, batch_size=10)
  list_names={}
  for method_name in methods.keys():
    add_to_list=[]
    for i in range(1,n_components+1):
      add_to_list.extend([method_name+feature_name+str(i)])
    list_names[method_name]=add_to_list

  all_results=pd.DataFrame()
  for i, (label, method) in enumerate(methods.items()):
    print("NOW RUNNING: ",label)
    Y = method.fit_transform(df)
    print(Y)
    df=pd.DataFrame(data=Y, index=None, columns=list_names[label]) 
    all_results=pd.concat([all_results,df], axis=1)
  with open(os.path.join(run_dir,"info.txt"), 'a') as file:
      file.write(""" The {}-encoded vector is reduced in dimensionality using {}\n""".format(feature_name, methods.keys()))

  return(all_results)

def PCA_reduction(X, nr_features):
  ipca = IncrementalPCA(n_components=nr_features, batch_size=10)
  new_X = ipca.fit_transform(X)
  return new_X



#------------------------------------------------------------------------------------------------------------------------------------
#FEATURE SELECTION
def feature_selection(X,y, RMSD_mode,nr_features, run_dir,manual_sele=None):
  feature_columns=X.columns
  all_target_cols=y.columns
  target_cols=[col for col in y.columns if "nom" in col]
  print("target_cols", target_cols)

  X_train, X_test, y_train, y_test = train_test_split(X, y[target_cols], test_size=0.33, random_state=42)
  X_train.reset_index(drop=True, inplace=True)
  y_train.reset_index(drop=True, inplace=True)
  #print(y_train.columns)
  for i in y_train.columns:
    if 'nom' in i:
      y_train=y_train[i]
      break
  print(y_train, X_train)

  selection_outcomes=collections.OrderedDict()
  feature_df=X
  feature_df_inuse=feature_df

  if manual_sele!=None:
    if any(isinstance(l, list) for l in manual_sele):
      for n , features in zip(range(len(manual_sele)), manual_sele):
        selection_outcomes["manual_"+str(n)]=features
    else:
      selection_outcomes["manual"]=manual_sele
   
  selected_features=random_forest_features(X_train, y_train, feature_df_inuse)
  selection_outcomes["random_forest_features"]=selected_features
  feature_df_inuse=feature_df
  print("random forest feature selection done")

  selected_features=anova(X_train,y_train, nr_features)
  selection_outcomes["anova"]=selected_features
  feature_df_inuse=feature_df
  print("anova done")

  selected_features=variance_threshold(feature_df_inuse)
  selection_outcomes["VarianceThreshold"]=selected_features
  feature_df_inuse=feature_df
  
  with open(os.path.join(run_dir,"info.txt"), 'a') as file:
    file.write("Feature selection outcomes for each method: {}\n".format(selection_outcomes))

  """selected_features=recursive_elimination(X_train,y_train, feature_df_inuse,run_dir=run_dir)
  selection_outcomes["recursive_elimination"]=selected_features
  print("recursive elimination done")"""

  commonly_selected_features=list(set.intersection(*[set(i) for i in list(selection_outcomes.values())]))
  selection_outcomes["commonly_selected"]=commonly_selected_features

  with open(os.path.join(run_dir,"info.txt"), 'a') as file:
    file.write("Commonly selected features: {}\n".format(commonly_selected_features))
    file.write("Different selection method outputs: {}\n".format(selection_outcomes))
  print(selection_outcomes)
  return selection_outcomes

def variance_threshold(X):
  columns = X.columns.values.tolist()
  thresh=(.8 * (1 - .8))
  selector = VarianceThreshold(threshold=thresh)  #change threshold value, is usually 0
  X=X.to_numpy()
  new_X=selector.fit_transform(X)
  labels = [columns[X] for X in selector.get_support(indices=True)]
  return labels
  
def calc_vif_score(X,del_y):
  feature_columns=([col for col in X.columns if col not in del_y])
  X=X[feature_columns]
  vif = pd.DataFrame()
  vif['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
  vif['variable'] = X.columns
  return vif

def random_forest_features(X_train, y_train, X):
  # fitting the model
  y_train=y_train.astype(str)
  sel = SelectFromModel(RandomForestClassifier(n_estimators=100, random_state=42, n_jobs=1))
  sel.fit_transform(X_train,y=y_train)
  sel.get_support()
  selected_features= X_train.columns[(sel.get_support())]
  return selected_features

def anova(X_train,y_train, nr_features):
  selector = SelectKBest(k=10)  #have to optimize the number of features to keep (usually 10)
  y_train=y_train.astype(str)
  new_X = selector.fit_transform(X_train, y=y_train)
  selector.get_support()
  selected_features= X_train.columns[(selector.get_support())]
  return selected_features

def recursive_elimination(X_train,y_train, feature_df_inuse,run_dir=None):
  if not os.path.exists(os.path.join(run_dir,"graphs")):
    os.makedirs(os.path.join(run_dir,"graphs"))
  visualizer_out_path_sub=os.path.join(run_dir,"graphs")
  features=X_train.columns
  X_normalized = preprocessing.normalize(X_train, norm='l2')
  y_train=y_train.astype(str)
  estimator = SVC(kernel="linear")
  mcc_scorer=make_scorer(matthews_corrcoef)
  estimator=estimator.fit(X_train,y=y_train)
  kfold = KFold(n_splits=10, shuffle=True)
  selector = RFECV(estimator, scoring=mcc_scorer,step=2, cv=kfold)
  selector.fit(X_train, y=y_train)
  selector.finalize()
  plt.savefig(os.path.join(visualizer_out_path_sub,"recursive_elimination.png"))
  plt.clf()

  ranked_features=selector.ranking_
  print(ranked_features)
  indexs= np.where(ranked_features == 1)
  selected_features=features[indexs]
  print(selected_features)
  return selected_features


if __name__ == '__main__':
  file=pd.read_csv(sys.argv[1])
  aa_encoding_all(file, ['tip_res', 'tip', 'res', 'atom_name', 'chain', 'pos', 'critical_seq'])
  #blosum_pam_score("H", "H", "BLOSUM35")