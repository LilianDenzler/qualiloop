#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Library
   File:       myfunctions.py
   
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
import os
import pandas as pd
import numpy as np
import joblib

from qualiloop import save_RMS_lib
from qualiloop import get_abYmod2_lib
from qualiloop import seq_parser_lib
from qualiloop import charge_lib
from qualiloop import happiness_lib
from qualiloop import accessibility_lib
from qualiloop import protrusion_lib
from qualiloop import get_voids_lib
from qualiloop import ecalc_lib
from qualiloop import critical_res_lib
from qualiloop import res_in_range_lib
from qualiloop import blosum_pairer_lib

import doctest




def get_abymod_input(upload_dir,ID, file):
	"""# return csv of input sequence for abymod
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> get_abymod_input(upload_dir, input_seq_df)#doctest: +NORMALIZE_WHITESPACE 
	'/home/lilian/sync_project/CDRH3lib/tests/abymod_input.csv'
	"""
	input_seq_df=get_input_seq(file)
	input_seq_df["position"] = input_seq_df["chain"] + input_seq_df["pos"]
	abymod_input_df=input_seq_df[['position', 'res']].copy()
	outfile = str(os.path.join(upload_dir, ID+"_input.csv"))
	abymod_input_df.to_csv(outfile, index = False, header = False, sep = '	', encoding = 'utf-8')
	return (str(outfile))



def get_models (upload_dir, ID, file, Redundant_PDBs=None):
	"""# return csv of input sequence for abymod
	>>> current="/home/lilian/sync_project/CDRH3lib/tests"
	>>> upload_dir="/home/lilian/sync_project/CDRH3lib/tests"
	>>> ID ="test_model"
	>>> test_file=os.path.join(current,"1F11_1.model")
	>>> input_seq_df=get_input_seq(test_file)
	>>> abymod_input_path= get_abymod_input(upload_dir, input_seq_df)
	>>> get_models(abymod_input_path,upload_dir, ID)
	('/home/lilian/sync_project/CDRH3lib/tests/test_model.pdb.model', '/home/lilian/sync_project/CDRH3lib/tests/test_model.log')
	"""
	abymod_input_csv=get_abymod_input(upload_dir, ID,file)
	print(abymod_input_csv)
	(model_path,log_path)=get_abYmod2_lib.pass_commands(abymod_input_csv, upload_dir, ID,Redundant_PDBs=Redundant_PDBs)
	return(model_path,log_path)


###########################################################################################################################


###########################################################################################################################


"""************
	This function takes the user input PDB model file and creates a dataframe of the sequence data. 

	###change to include more of the PDB values
	***************"""
def get_input_seq(file):
	file=str(file)
	input_seq_df=seq_parser_lib.get_df(file)
	return (input_seq_df)

"""************
	Outputs df of sequence of input PDB model file. 
	***************"""
def get_loop_seq (file):
	input_seq_df=get_input_seq(file)
	columns=["pos","res"]
	copy = False
	seq=[]
	pos=[]
	for index, row in input_seq_df.iterrows():
		if row["chain"]=="H" and row["pos"]=="95":
			copy = True
			seq+=row["res"] 
			pos+=[row["pos"]]
			continue
		elif row["chain"]=="H" and row["pos"]=="102":
			seq+=row["res"] 
			pos+=[row["pos"]]
			copy = False
			continue
		elif copy:
			seq+=row["res"] 
			pos+=[row["pos"]]
	loop_seq_df=pd.DataFrame({'pos':pos, 'res':seq},columns=columns, index=None)
	length=loop_seq_df.shape[0]
	seq_string = loop_seq_df['res'].tolist()
	seq_string="".join(seq_string)
	length_seq_df=pd.DataFrame({'seq':seq_string, 'length':[length]},columns=["seq","length"], index=None)
	return (loop_seq_df,length_seq_df)


"""************
	Outputs df of loop total charge, nr of charged residues in the loop.
	***************"""

def get_loop_charge(file):
	loop_seq_df,length_seq_df=get_loop_seq (file)
	charge_df=charge_lib.get_charge(loop_seq_df)
	return (charge_df)


"""************
	Outputs df of Happiness_score
	***************"""
def get_happiness_score(file):
	happiness_df=happiness_lib.get_happy(file)
	return (happiness_df)


"""************
	Outputs df of loop protrusion (i.e distance in angstrom from line between base residues and tp residue.)
	***************"""
def get_protrusion(file):
	protrusion_df=protrusion_lib.calc_protrusion(file)
	#protrusion_df["tip"] = protrusion_df["tip_res"] + protrusion_df["tip_pos"]
	#protrusion_df.drop(labels=["tip_res", "tip_pos"], axis=1)
	return (protrusion_df)



def similarity_identity(file, template=None):
	#not finished
	try:
		if template is not None:
			columns=["identity","similarity", "template"]
			data=identity,similarity,template
			df = pd.DataFrame(data, columns=columns)
			return (df)
	except:
		return pd.DataFrame()

"""************
	Outputs df of identity and similarity to template sequence. 
	***************"""
def handle_log_file(log_name):
	columns=["identity","similarity", "template", "target"]
	data=[]
	with open(os.path.join(log_name), "r") as infile:
		for line in infile:
			#INFO: CDR-H3 (YEIR/YEWA) SeqID: 0.500 Similarity: 0.381
			if "INFO: CDR-H3" in line:
				output=line.split(" ")
				template_target=output[2]
				template_target=template_target.split("/")
				target=template_target[0]
				template= template_target[1]
				target=target.replace("(", "")
				template=template.replace(")", "")

				identity=output[4]
				similarity=output[6]
				similarity=similarity.replace('"', "")
				similarity=similarity.replace("\n","")
				write=[identity,similarity, template, target]
				data.append(write)

	df = pd.DataFrame(data, columns=columns)
	return (df)




"""************
	Outputs df of Hydrophobicity values based on the consensus values by Eisenberg et al. 
	(Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109). 
	output:
	Mean of hydrophobicity values of loop
	Sum of absolute Differences between loop and template loop
	***************"""				
def hydropathy(file,template=None):
	(loop_seq_df,length_seq_df) =get_loop_seq(file)
	#Consensus values: Eisenberg, et al 'Faraday Symp.Chem.Soc'17(1982)109
	Hydrophathy_index = {'A': 00.250, 'R': -1.800, "N": -0.640, "D": -0.720, "C": 00.040, "Q": -0.690, "E": -0.620, "G": 00.160, "H": -0.400, "I": 00.730, "L": 00.530, "K": -1.100, "M": 00.260, "F": 00.610, "P": -0.070,
							"S": -0.260, "T": -0.180, "W": 00.370, "Y": 00.020, "V": 00.540, "X": -0.5}#-0.5 is average
	
	hydro_value=0
	data=[]
	for index, row in loop_seq_df.iterrows():
		aminoacid=row['res']
		hydro_value+=Hydrophathy_index[aminoacid]

	if template==None:
		columns=["Hydropathy"]
		write=[hydro_value]
	else:
		columns=["Hydropathy","Hydropathy_diff"] 
		target=length_seq_df["seq"].iloc[0]
		for a,b in zip(template, target):
			hydro_diff=abs(abs(Hydrophathy_index[b])-abs(Hydrophathy_index[a]))
		write=[[hydro_value,hydro_diff]]

	df = pd.DataFrame(write, columns=columns, index=None)

	return (df)


"""************
	Outputs df of accessibility 
	***************"""				   
def accessibility(file):
	df=accessibility_lib.pdbsolv_run(file)
	return (df)

def bulged_non_bulged(feature_directory,input_seq_directory,model_directory,bioptools_directory):
	pass

def seq_ransomisation(feature_directory,input_seq_directory):
	pass
	#also information entropy

def voids(file):
	voids_df=get_voids_lib.get_voids(file)
	return voids_df
	
def get_in_range(file):
	(loop_seq_df,length_seq_df)=get_loop_seq(file)
	in_range_df=res_in_range_lib.get_in_range(file)
	length=length_seq_df["length"].astype('float').to_numpy()[0]
	contacts_all_norm=in_range_df["contacts_all"].astype('float').to_numpy()[0]/length
	contacts_out_norm=in_range_df["contacts_out"].astype('float').to_numpy()[0]/length
	in_range_df["contacts_all_norm"]=contacts_all_norm
	in_range_df["contacts_out_norm"]=contacts_out_norm
	return in_range_df
	

def get_ecalc(file):
	energy_df=ecalc_lib.get_ecalc(file)
	return energy_df

def gromacs_energy(log_name):
	energy_df = pd.DataFrame()
	energy=None
	with open(os.path.join(log_name), "r") as infile:
		for line in infile:
			#Potential Energy  = -5.2766569e+05
			if "Potential Energy" in line:
				output=line.split("= ")[1]
				energy=float(output)
	if energy==None:
		return pd.DataFrame()
	energy_df["energy"]=[energy]
	return energy_df

	


def critical_res(file):
	(loop_seq_df,length_seq_df)=get_loop_seq(file)
	critical_df=critical_res_lib.get_critical_res(loop_seq_df)
	return critical_df

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

def similarity_length(file, template=None, log_path=None):
	try:
		if log_path is not None:
			log_df=handle_log_file(log_path)  # log_df outputs columns: "identity","similarity", "template", "target"
			template=log_df["template"].iloc[0]
			similarity_df=log_df[["identity","similarity"]]
			similarity=similarity_df["similarity"].astype('float').to_numpy()[0]
		elif template is not None:
			similarity_df=similarity_identity(file, template)
			similarity=similarity_df["similarity"].astype('float').to_numpy()[0]
		else: 
			similarity_df=None
			return None
	except:
		similarity_df=None
		return None

	loop_seq_df,length_seq_df=get_loop_seq(file)
	length=length_seq_df["length"].astype('float').to_numpy()[0]
	simlength=similarity/length
	simlength_df=pd.DataFrame()
	simlength_df["simlength"]=[simlength]
	print(simlength_df)
	return simlength_df


def blosum_pair(file,template=None):
	(loop_seq_df,length_seq_df)=get_loop_seq(file)
	result_df=blosum_pairer_lib.run_pair_scorer(file,loop_seq_df,template=template)
	return result_df



def merge(file,list_of_features=None, template=None, log_path=None):
	final_df=[]
	model_features=[]
	#--------------------------------------------------------------------------------
	#calculate all the separate feature dataframes
	#--------------------------------------------------------------------------------
	# get identity and similarity from either the log-file or the template sequence
	try:
		if log_path is not None:
			log_df=handle_log_file(log_path)  # log_df outputs columns: "identity","similarity", "template", "target"
			template=log_df["template"].iloc[0]
			similarity_df=log_df[["identity","similarity"]]
		elif template is not None:
			similarity_df=similarity_identity(file, template)
		else: 
			similarity_df=None
	except:
		template=None
		similarity=None
		similarity_df=pd.DataFrame()
		log_df=pd.DataFrame()
	#get length
	(loop_seq_df,length_seq_df)=get_loop_seq(file)
	#get simlength
	simlength_df=similarity_length(file, template=template, log_path=log_path)
	#get total charge and nr. charged, 
	charge_df=get_loop_charge(file)
	#get happiness score
	happiness_df=get_happiness_score(file) 
	        
	#get protrusion in Angstrom and tip residue residue and position 
	protrusion_df=get_protrusion(file)
	
	#get total hydropathy if not template, also get difference in hydropathy if template is provided
	try:
		if template is not None:
			hydropathy_df=hydropathy(file,template)
		else: 
			hydropathy_df=hydropathy(file)
	except:
		hydropathy_df=pd.DataFrame()
	#get accessibility values (take averages of residues and totals of: 
	#accessibility, relative accessibility, side-chain accessibility, relative side-chain accessibility)
	
	accessibility_df=accessibility(file)
	
	#	accessibility_df=pd.DataFrame()
	#get voids and nearest atom coords
	"""try:
		voids_df=voids(file)
	except:
		voids_df=pd.DataFrame()"""
	#critical residues
	try:
		critical_df=critical_res(file)
	except:
		critical_df=pd.DataFrame()
	
	try:
		energy_df=gromacs_energy(log_path)
	except:
		energy_df=pd.DataFrame()
	
	#get energies
	#energy_df=get_ecalc(file) #########################################################################''''''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	blosum_pair_df=blosum_pair(file,template=template)
	#except:
	#	blosum_pair_df=pd.DataFrame()

	#try:
	in_range_df=get_in_range(file)
	#except:
	#	in_range_df=pd.DataFrame()

	ID_df=pd.DataFrame()
	ID_code=str(file.split(".")[0].split("/")[-1])
	ID_df["ID"]=[ID_code]


	#--------------------------------------------------------------------------------------------
	#concatenate all the separate df into one single-row dataframe
	#--------------------------------------------------------------------------------------------
	append_list=[similarity_df, length_seq_df, simlength_df, charge_df, happiness_df, protrusion_df, 
	hydropathy_df, accessibility_df, critical_df, blosum_pair_df, in_range_df, energy_df, ID_df] ###voids takes toooo loooong

	for i in append_list:
		print (i)
	model_features = pd.concat(append_list, axis=1)
	if list_of_features==None:
		return model_features
	if list_of_features:
		model_features=model_features[list_of_features]
		#sanity check:
		if model_features.columns!=list_of_features:
			raise Exception("The list of list of features is formatted incorrectly or contains features that are not included in Qualiloop") 
		else:
			return model_features

def run(file,list_of_features=None,template=None, log_path=None):
	model_features= merge(file,list_of_features, template, log_path)
	return model_features
	



if __name__ == '__main__':
	upload_dir=sys.argv[1]
	Redundant_PDBs=sys.argv[3]
	for file in os.listdir(sys.argv[2]):
		ID=file.replace(".pdb", "")
		print(ID)
		file=os.path.join(str(sys.argv[2])+file)
		model_path,log_path=get_models(upload_dir, ID, file,Redundant_PDBs=Redundant_PDBs)
	

	"""file=str(os.path.join("/home/lilian/sync_project/CDRH3lib/tests","1F11_1.pdb"))
	log_path=os.path.join("/home/lilian/sync_project/CDRH3lib/tests","1F11_1.log")
	result=run(file, log_path=log_path)"""
	