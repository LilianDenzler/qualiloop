import os
import pandas as pd
import numpy as np
import sys
from sklearn.decomposition import TruncatedSVD
import collections
from sklearn import manifold, datasets
from functools import partial


def get_critical_res(loop_seq_df):
    #Critical 6 residues (H94-H97, H101-H102)
    #critical_seq1=(loop_seq_df.loc[(loop_seq_df['pos'] == "94"),"res"]).to_list()
    critical_seq2=(loop_seq_df.loc[(loop_seq_df['pos'] == "95"),"res"]).to_list()
    critical_seq3=(loop_seq_df.loc[(loop_seq_df['pos'] == "96"),"res"]).to_list()
    critical_seq4=(loop_seq_df.loc[(loop_seq_df['pos'] == "97"),"res"]).to_list()
    critical_seq5=(loop_seq_df.loc[(loop_seq_df['pos'] == "101"),"res"]).to_list()
    critical_seq6=(loop_seq_df.loc[(loop_seq_df['pos'] == "102"),"res"]).to_list()
    critical_seq=""
    try:
        for i in [critical_seq2,critical_seq3, critical_seq4, critical_seq5, critical_seq6]:
            critical_seq+=i[0]
        print("SEQ", critical_seq)
        print("SEQ", critical_seq)
        #Critical 6 residues (H94-H97, H101-H102)
        critical_df=pd.DataFrame()
        critical_df["critical_seq"]=[critical_seq]
        return (critical_df)
    except:
        return pd.DataFrame()

    #python3 myfunctions.py ~/sync_project/input_Abymod_seq/

    #6x20 Integer (1-hot)6x20 
    #Real (BLOSUM80)6x20 
    #Real (BLOSUM62)6x20 
    #Real (BLOSUM45)6x20 
    #Real (PAM250/MDM)6x5 
    #Real (Li-Koehl PCA5)6x4 
    #Real (Li-Koehl PCA4)6x3 
    #Real (Li-Koehl PCA3)6x3 
    #Real (El-Maarty 3-parameter physical encoding)6x4 
    #Real (Abhinandan 4-parameter physical encoding)
