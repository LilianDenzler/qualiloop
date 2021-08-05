#!/usr/bin/python3
"""/*************************************************************************
   Program:    Blosum-pair scorer graph maker
   File:       blosum_pairer_graph.py
   
   Version:    V.1
   Date:       19.01.21
   Function:   General type definitions, defines, macros and globals
   
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
   The script is used to generate a feature for the CDRH3 loop model quality 
   predictor. 
   1. take first res n at one end and pair with each n+2,n+3,...until the other 
   end of loop
   2. for each pair calculate the sum of the blosum score difference between the
   model's residue and the template's residue at this position. 
   3. Take the lowest scoring)i.e. most neagtive residue pair
   4. calculate the -log2 of the number of residue between the lowest scoring
   pair
   5. add the lowest score (which is neg) to the -log2(nr. of res distance)
   6. look for correlation between the value from 5 and the RMSD of the model

   If there is correlation-> make it a feature
   maybe adjust weighting between distance term and loswest blosum score term. 
**************************************************************************
   Usage:
   ======
   This library is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   20.11.20 Original 
 
*************************************************************************/"""
import sys
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./CDRH3lib')
import myfunctions
import preprocessor
import modelmaker
import sys
import os
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
import blosum_pairer_lib

def full_run(model_dir, actual_dir):
   final_df=pd.DataFrame()
   for file in os.listdir(model_dir):
      model_path=os.path.join(model_dir,file)
      filename=file.replace(".pdb.model", "")
      actual_path=os.path.join(actual_dir, filename+".pdb")
      #print(model_path, log_path, actual_path)
      if not os.path.isfile(model_path) or not os.path.isfile(actual_path):
         continue

      result_df=myfunctions.get_in_range(model_path)
      print(result_df)
      
      RMSD_num_df=modelmaker.get_RMSD_num(model_path, actual_path)
      #print(result_df)
      row=pd.concat([RMSD_num_df, result_df], axis=1)
      final_df=pd.concat([final_df, row])
      #input(final_df)
      
   print(final_df)
   final_df.to_csv(os.path.join("./in_range.csv"))
   return final_df

def plot(final_df):
   for y_name in ("local_AA", "local_CA", "global_AA", "global_CA"):
      for x_name in ("contacts_all", "contacts_out", "contacts_ratio_out_all"):
         data=final_df[[x_name, y_name]]
         sns.jointplot(data=data, x=x_name, y=y_name, kind='reg',stat_func=stats.pearsonr,joint_kws={'line_kws':{'color':'black'}},scatter_kws={'alpha':0.0})
         ax=sns.scatterplot(data=data, x=x_name, y=y_name,alpha=0.5, s=50)
         plt.savefig(os.path.join("./in_range_graphs/","in_range"+x_name+y_name+".png"))
         plt.show()


if __name__ == '__main__':
   #final_df=full_run(sys.argv[1], sys.argv[2])
   path=sys.argv[1]
   print(path)
   final_df=pd.read_csv(path, header=0)
   final_df=final_df[["contacts_all", "contacts_out", "contacts_ratio_out_all", "local_AA", "local_CA", "global_AA", "global_CA"]]
   print(final_df)
   plot(final_df)