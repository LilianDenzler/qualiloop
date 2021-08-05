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
  This program replicates the technique for amino-acid encoding used in the 
  paper by Martin & Abhinandan 2010. The amino acids were encoded in 4d vectors 
  for neural network input preparation.The four features calculated were:
  #1. total number of side-chain atoms
  #2. number of side-chain atoms in shortest path from Calpha to most distal atom
  #3. eisenberg consensus hydrophobicity 
  #4. charge (histidine was assigned +0.5)
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
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./CDRH3lib')
sys.path.append('~/sync_project/WWW/CDRH3loop')
import os
import pandas as pd
import numpy as np
import joblib
import time

def run(seq):
	columns=["res_charge", "res_sc_nr", "res_compactness", "res_hydrophob"]
	seq_df=pd.DataFrame(columns=columns)
	for res in seq:
		charge_of_res=charge(res)
		hydrophobicity_res=hydophobicity(res)
		compactness_res=compactness(res)
		nr_side_chain_atoms_res=nr_side_chain_atoms(res)

		seq_df=seq_df.append({"res_charge":charge_of_res, "res_sc_nr":nr_side_chain_atoms_res, "res_compactness":compactness_res, "res_hydrophob": hydrophobicity_res}, ignore_index=True)

	return seq_df

def nr_side_chain_atoms(res):
	#1. total number of side-chain atoms
	nr_side_chain_atoms_dic = {'A': 1, 'R': 7, "N": 4, "D": 4, "C": 2, "Q": 5, "E": 5, "G": 0, "H": 6, "I": 4, 
								"L": 4, "K": 15, "M": 4, "F": 7, "P": 4,
								"S": 2, "T": 3, "W": 10, "Y": 8, "V": 3 } #"X": 10.375
	nr_side_chain_atoms=nr_side_chain_atoms_dic[res]
	return nr_side_chain_atoms
	


def compactness(res):
	#2. number of side-chain atoms in shortest path from Calpha to most distal atom
	compactness_dic = {'A': 1, 'R': 6, "N": 3, "D": 3, "C": 2, "Q": 4, "E": 4, "G": 0, "H": 4, "I": 3, 
						"L": 3, "K": 6, "M": 4, "F": 5, "P": 2,
						"S": 2, "T": 2, "W": 6, "Y": 6, "V": 2} #, "X": 4.45
	compactness=compactness_dic[res]
	return compactness
	
	

def hydophobicity(res):
	#3. eisenberg consensus hydrophobicity 
	#Consensus values: Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109
	Hydrophathy_index = {'A': 00.250, 'R': -1.800, "N": -0.640, "D": -0.720, "C": 00.040, "Q": -0.690, "E": -0.620, "G": 00.160, "H": -0.400, "I": 00.730, "L": 00.530, "K": -1.100, "M": 00.260, "F": 00.610, "P": -0.070,
							"S": -0.260, "T": -0.180, "W": 00.370, "Y": 00.020, "V": 00.540, "X": -0.5}#-0.5 is average
	
	hydrophobicity=Hydrophathy_index[res]
	return hydrophobicity
	
def charge(res):
	dic = {"D":-1, "K": 1,"R": 1,'E': -1, 'H':0.5}
	charge=0
	if res in dic:
		charge+=dic[res]
	return (charge)
