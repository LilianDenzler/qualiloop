#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Deep Feature Generator
   File:       deep_feature_gen.py
   
   Version:    V1.1
   Date:       22.03.21
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
   The program contains all functions needed for generating new features from
   an existing feature datafram in a process called deep feature generation.


**************************************************************************
   Usage:
   ======
   This library is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   22.03.21 Original
*************************************************************************/"""
import sys
import os
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./')
sys.path.append('~/sync_project/WWW/CDRH3loop')
import numpy as np
import pandas as pd
import featuretools as ft


def feature_synthesis():
  #do this for blosum63 dfs and NFL and other matrices of all loops!
  #also do for finished dataset
  #ft.list_primitives() run to see all features
  data=pd.read_csv("~/sync_project/Feature/features+critical2.csv")
  target_df=data[["local_CA", "ID"]]
  del data["local_AA"]
  del data["local_CA"]
  del data["global_AA"]
  del data["global_CA"]
  #data['ID'] = "first_df"
  es = ft.EntitySet(id = 'data')
  es.entity_from_dataframe(entity_id = 'data', dataframe = data, index="ID")
  print(es)

  es.entity_from_dataframe(entity_id="local_CA",dataframe=target_df, index="ID")
  print(es)

  feature_matrix, feature_defs = ft.dfs(entityset=es,target_entity="data",max_depth = 2, verbose = 1, n_jobs = -1)
  #es.entity_from_dataframe(entity_id="data",dataframe=target_df,index="ID")
  #feature_matrix, feature_names = ft.dfs(entityset=es,target_entity = 'local_CA', max_depth = 2, verbose = 1, n_jobs = 3)
  return feature_matrix, feature_defs


if __name__ == '__main__':
  possible=ft.list_primitives()
  possible=possible
  feature_matrix, feature_defs=feature_synthesis()
  print(feature_matrix)
  print(feature_defs)