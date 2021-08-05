#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Bulged Loop Finder
   File:       bulged_lib.py
   
   Version:    V2.1
   Date:       07.01.21
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
   The library contains all functions needed to extract features of the input 
   PDB model file. A dataframe containing all the extracted features will be outputted.
**************************************************************************
   Usage:
   ======
   This library is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   20.11.20 Original - library for making df of full dataset for testing, 
			inclusing simple features such as charge, seq ID, similarity, length
   V0.1   15.12.20 Update  -  library more modular, incorporating helper.py files
			instead of all functions in one script
   V1.1	  02.02.21 Update  -  library incorporates more features

   V2.1   07.01.21 Update  -  library incorporates more features, more modular
			restructured, so so it can be used for model training 
			(i.e. use full dataset) and for running the model with unknown (i.e. user)

*************************************************************************/"""
import sys
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./CDRH3lib')
import os
import pandas as pd
import numpy as np
import joblib



'''Previous studies identified a sequence motif (Arg or Lys at T2 and Asp at T6) that contributes to bulged torso formation in some but not all cases;
	these key residues were conserved in our bulged cluster, with 80% of bulged structures presenting Arg or Lys at T2, 73% presenting Asp at T6, and 65%
	retaining the complete T2/T6 sequence motif (S1 Fig) [9,10]. - Finn et al. PLOS One, 2016'''

'''C-terminal loop residues form a pseudo dihedral angle,a101(100X–103 using Chothianumbering),
	of 39 ̊ and a pseudo bond angle,t101(100X–102using Chothia numbering), of 101 ̊

	They go on to use 3 standard deviations, so ca: pseudo dihedral a101=>6.5-70.7° and pseudo bond angle t101 =>81.81-117.81°
	t101, the Ca–Ca–Ca pseudo bond angle for the three C-terminal residues;and 2)a101,the Ca–Ca–Ca–Ca pseudo dihedral angle for
	the three C-terminal residues in the CDR H3loop and one adjacent residue inthe H chain framework

	-Weizner and Gray, 2016, Journal of Immunology'''

'''for filename in os.listdir(input_seq_directory):
		file=open(os.path.join(input_seq_directory,filename))
		for line in file:
			line.split(" ")
			if line[1]=="86":'''
pwd=os.getcwd()
os.chdir(bioptools_directory)
for filename in os.listdir(model_directory):
	name= filename.replace(".pdb.model","")
	angles=os.popen("pdbtorsions -c {}".format(os.path.join(model_directory,name+".pdb"+".model"))).readlines()
	for line in angles:
		columns=line.split()
		print (line)
		if "H100" in columns[0]:
			angle100x=columns[2]
			AAangle100x=columns[1]
			if float(angle100x) >=6.5 and float(angle100x) <= 70.7:
				angle100x="k"
			else:
				angle100x="e"
		if "H101" in columns[0]:
			angle101=columns[2]
			AAangle101=columns[1]
			if float(angle101) >=6.5 and float(angle101) <= 70.7:
				angle101="k"
			else:
				angle101="e"
		if "H102" in columns[0]:
			angle102=columns[2]
			AAangle102=columns[1]
			if float(angle102) >=6.5 and float(angle102) <= 70.7:
				angle102="k"
			else:
				angle102="e"
		if "H103" in columns[0]:
			angle103=columns[2]
			AAangle103=columns[1]
			if float(angle103) >=6.5 and float(angle103) <= 70.7:
				angle103="k"
			else:
				angle103="e"
			write=[angle100x, angle101, angle102, angle103, name]
			columnids=["angle100x","angle101","angle102","angle103","ID"]
			if os.path.exists(os.path.join(feature_directory,"kinked"+ ".csv")): #use feature_directory,"RMS_by_res_feature_csv_"+dics2[y] + ".csv") for separate files
				with open(os.path.join(feature_directory,"kinked" + ".csv"), "a") as f:
					writer = csv.writer(f)
					writer.writerow(write)
			else:
				with open(os.path.join(feature_directory,"kinked"+ ".csv"),"w") as f:
					writer = csv.writer(f)
					writer.writerow(columnids)
					writer.writerow(write)