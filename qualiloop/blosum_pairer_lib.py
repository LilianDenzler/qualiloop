#!/usr/bin/python3
"""/*************************************************************************
	 Program:    Blosum-pair scorer
	 File:       blosum_pairer_lib.py
	 
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
import sys
import os
import pandas as pd
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from scipy import stats
import glob

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


	#print(matrices_path)
	return matrices_path[0]

def blosum_pam_score(res1, res2, matrix):
	matrices_path=get_file_paths(matrix)
	#filepath=os.path.join((str(os.path.abspath("data/matrices_blosum_pam"))), matrix)
	matrix=[]
	with open(matrices_path) as file:
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

def pair_scoring(file, matrix,loop_seq_df,template=None):
	if template==None:
		return pd.DataFrame()
	score_list=[]
	pair_scores=[]
	pos_list=[]
	counter=0
	loop_positions=loop_seq_df["pos"].tolist()
	if len(template)!= len(loop_positions) or len(loop_positions)==0:
		return pd.DataFrame()
	for count in range(0, len(template)):
		temp_res=template[count]
		row=loop_seq_df.iloc[count]
		pos=row["pos"]
		res=row["res"]
		score_N=blosum_pam_score(temp_res, res, matrix)
		score_list.append(score_N)
		pos_list.append(pos)
	for a in range(0, len(pos_list)):
		for b in range(a+2, len(pos_list)):
			pair_score=int(score_list[a])+int(score_list[b])
			res_a=pos_list[a]
			res_b=pos_list[b]
			pair_scores.append([pair_score,res_a, res_b])
	 
	worst_score=min([x[0] for x in pair_scores])
	 
	pair_scores_df=pd.DataFrame(pair_scores, columns=["score", "res_a", "res_b"])
	res_a_worst=pair_scores_df.loc[pair_scores_df['score'] == worst_score][["res_a"]].to_numpy()[0]
	res_b_worst=pair_scores_df.loc[pair_scores_df['score'] == worst_score][["res_b"]].to_numpy()[0]
	index_a=loop_positions.index(res_a_worst)
	index_b=loop_positions.index(res_b_worst)
	distance=abs(index_a-index_b)
	distance_log=(math.log(distance,2))*(-1)
	result=distance_log*worst_score
	result_df=pd.DataFrame()
	result_df["blosum_dist"]=[result]
	return result_df
	 
def run_pair_scorer(file, loop_seq_df,template=None):
	final_df= pair_scoring(file, "BLOSUM80",loop_seq_df,template=template)

	return final_df
