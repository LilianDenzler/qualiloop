#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Model Maker
   File:       modelmaker.py
   
   Version:    V1.1
   Date:       09.03.21
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
   The program contains all functions needed to calculate various target y 
   values that the model will predict. It also builds various different models 
   and assesses their performance. The types of models that will be assessed 
   can be restricted by making specifications in the config file. If no such 
   specifications made, all possible models in this program will be assessed 
   and the most succesful one will be saved. 

   The config file takes the following input:
   - path to save the most optimal model
   - simple or structured model (see documentation for further information)
	 (if not specified, both simple and structures models will be assessed)
   - if a log file is inputted or not by the user using the final model 
	 (default is yes, there is a log file)

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
import math 

from qualiloop import myfunctions
from qualiloop import save_RMS_lib

#import staggered_bins_lib #?
from qualiloop.hyperparam_opti import toplayer as opti_gridsearch
from qualiloop.hyperparam_opti_genetic import opti_run as opti_genetic
#import visualizer

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
from sklearn.metrics import make_scorer
import xgboost as xgb

import doctest
import gc

#from pyinstrument import Profiler
#python3.6 modelmaker.py ~/sync_project/abymod_structures2 ~/sync_project/abymod_structures2_log/ ~/sync_project/actual_PDBs_NR/ /serv/www/html_lilian/RUDIMENTARY_DF.csv
#################################################################################################################################################################
#Functions for creating dataframe with which to train model


def get_RMSD_num(file, actual_file, mode):
	"""
	Calculate the numerical RMSD of the model and actual file in the loop 
	region using ProFit.
	The RMSD will be made for local (i.e. only loop aligned) and global 
	(i.e. whole structure , except loops aligned)and for only Ca-Ca and atom-atom each. 

	Input:  file      --- The model file in PDB format
			actual_file ---  The 'real' PDB file taken from the AbDb database

	Return: RMSD_num_df --- Dataframe containing the columns 
							"local_AA","local_CA","global_AA","global_CA"

	"""
	RMSD_num_df=save_RMS_lib.run_num(file,actual_file, mode)
	#print(RMSD_num_df,"RMSD_num_df")
	return RMSD_num_df


def get_all_bins(RMSD_mode, file,actual_file, threshold_list, RMSD_num_df=None):
	"""
	Calculate the binary RMSD threshold values for the file. 

	The RMSD_num_df argument was added to save time, as if get_RMSD_num was called before, 
	the second calculation of RMSD_num_df in this function can be avoided. 

	Input:  file      --- The model file in PDB format
			actual_file ---  The 'real' PDB file taken from the AbDb database
			threshold_list ---  list of threshold values in Angstroms. Is the 
								RMSD is above or exactely this value, the bin value 
								will be 1, if the RMSd is 0. 
			RMSD_num_df (optional) --- dataframe of numerical local_AA, local_CA, global_AA, global_CA RMSD values

	Return: df_all_bins --- Dataframe containing the columns: 
							"local_AA_bin"+str(threshold),"local_CA_bin"+str(threshold),"global_AA_bin"+str(threshold), "global_CA_bin"+str(threshold)
							for each threshold specified in threshold_list
							e.g. "local_AA_bin1, ...global_CA_bin8"

	"""
	df_all_bins=save_RMS_lib.get_all_bins(RMSD_mode, file,actual_file, threshold_list, RMSD_num_df)
	print("df_all_bins",df_all_bins)
	return df_all_bins

def RMSD_nom(RMSD_mode, threshlow,threshup, file, actual_file,counter=None, RMSD_num_df=None):
	"""
	Calculate the nominal RMSD values for the file. If the RMSD value is smaller than a value 
	from upper threshold list and also greater than the corresponding value in the lower threshold 
	list, the value will be the number of this nominal category.

	The RMSD_num_df argument was added to save time, as if get_RMSD_num was called before, 
	the second calculation of RMSD_num_df in this function can be avoided. 

	Input:  file      --- The model file in PDB format
			actual_file ---  The 'real' PDB file taken from the AbDb database
			threshlow ---  List of lower threshold values of the nominal categories
			threshup ---  List of upper threshold values of the nominal categories
			counter (optional) ---  If multiple nominal classes are to be tested, this counter 
									value indicates what number this nominal class is.
									It will be written in the df titles.

			RMSD_num_df (optional) --- Dataframe of numerical local_AA, local_CA, global_AA, global_CA RMSD values

	Return: df_nom --- Dataframe containing the values corresponding to the nominal RMSD class the structure is in.
						e.g. if counter=2 the columns would be:
						local_AA_nom2,...global_CA_nom2 


	"""
	print(threshlow,threshup)
	df_nom=save_RMS_lib.RMSD_nom(RMSD_mode, threshlow,threshup, file, actual_file,counter, RMSD_num_df)
	return df_nom


def make_df_nombin(RMSD_mode, file, actual_file,threshold_list, threshup, threshlow, RMSD_num_df=None):
	"""
	Make a dataframe containing all RMSD-bin class and RMSD-nominal class values for a file.

	Input:  file      --- The model file in PDB format
			actual_file ---  The 'real' PDB file taken from the AbDb database
			threshold_list ---  list of bin threshold values in Angstroms. Is the 
								RMSD is above or exactely this value, the bin value 
								will be 1, if the RMSd is 0. 
			threshlist --- list of tuples containing threshlow and threshup 
							(threshlow ---  List of lower threshold values of the nominal categories)
							(threshup ---  List of upper threshold values of the nominal categories)

			counter (optional) ---  If multiple nominal classes are to be tested, this counter 
									value indicates what number this nominal class is.
									It will be written in the df titles.

			RMSD_num_df (optional) --- Dataframe of numerical local_AA, local_CA, global_AA, global_CA RMSD values

	Return: df_nom --- Dataframe containing the values corresponding to the nominal RMSD class the structure is in.
						e.g. if counter=2 the columns would be:
						local_AA_nom2,...global_CA_nom2 
	"""
	df_all_bins=get_all_bins(RMSD_mode, file,actual_file, threshold_list, RMSD_num_df=RMSD_num_df)
	counter=0
	
	df_nom=RMSD_nom(RMSD_mode, threshlow,threshup, file, actual_file,counter, RMSD_num_df=RMSD_num_df)
	#combine to one dataframe
	all_RMSD=df_all_bins.merge( df_nom,on=["ID"])
	print("RMSD_num_df",df_all_bins, df_nom)
	return all_RMSD


def more_features(file, log_name=None, template=None):
	"""
	Make a dataframe containing all features that can be extracted and calculated using 
	the qualiloop library.The program is designed to readily incorporate any new features 
	when made available.
	A full datafram with all features is calculated and then used for feature selection in 
	later stages.
	If a log file of the abymod modelling process is not provided, some important features such 
	as similarity and identity will not be calculated compromising model performance. 
	The log file path is specified in the config file.
	Alternatively the similarity and identity can also be given. 

	Input:  file      --- The model file in PDB format
			log_name(optional)  --- the abymod log file 
			identity(optional)   ---   the sequence identity between the template sequence 
										and the target loop sequence
			similarity(optional)   ---   the sequence similarity between the template sequence 
										and the target loop sequence
			template(optional)   ---   template sequence as string used in the modellign process
			upload_dir(optional)  ---  If an upload directory path is specified, then the finished 
										feature dataframe will be saved here under the name:
										os.path.join(upload_dir,"features"+ID+".csv") 
										(the ID being the basename of the file without pdb extension)
 

	Return: model_features  --- Dataframe containing all features that can be calculated 
								using the qualiloop library, given if a log file, identity, 
								similarity, or the template sequence is provided: 
							

	"""
	filename=os.path.splitext(os.path.basename(file))[0]
	filename=filename.replace(".pdb", "")
	print(file, template, log_name)
	model_features= myfunctions.merge(file,list_of_features=None, template=template, log_path=log_name)
	return model_features


def get_bin_nom_dic_balanced(balanced_params,RMSD_num_df, mode):
	"""
	Balanced Bin mode:
	takes number of balanced bins to make and the stepsize.
	Stepsize gives the number of bins that will be combined to form a class. 
	If there are 10 bins and the stepsize is 2, there will be 5 classes. 


	#TO ADD IN FUTURE: HAVE MIN AND MAX
	"""
	

	RMSD_thresh_vals={}
	for RMSD_mode, balanced_params in zip(mode, balanced_params):
		[bins_in_class_RMSD, nr_of_balanced_bins]=balanced_params
		index=int(math.floor(RMSD_num_df.shape[0]/(nr_of_balanced_bins+1)))
		RMSD_array=RMSD_num_df[RMSD_mode].to_numpy()
		RMSD_array = RMSD_array.astype(np.float)
		sorted_target=np.sort(RMSD_array)
		values=[]
		for i in range(1, nr_of_balanced_bins+1):
			value = sorted_target[index * (i)]
			values.append(float(value))
		RMSD_thresh_vals[RMSD_mode+"_values"]=values

		threshold_list_nom=values[::bins_in_class_RMSD]
		threshup=[]
		if threshold_list_nom[0]!=0:
			new_thresh=[0]
			threshlow=np.append(new_thresh, threshold_list_nom)
			threshup=threshlow[1]
		threshup=np.append(threshup,threshold_list_nom[1:])
		threshup=np.append(threshup,[100])
		RMSD_thresh_vals[RMSD_mode+"_threshold_list"]=values
		RMSD_thresh_vals[RMSD_mode+"_threshup"]=threshup
		RMSD_thresh_vals[RMSD_mode+"_threshlow"]=threshlow
	
	return RMSD_thresh_vals

def get_bin_nom_dic_staggered(staggered_params,mode):
	"""
	Setting the binary and class thresholds when using the staggered mode. 

	Input:  staggered_params      ---   [staggering_value, max_angstrom, first_layer_size, second_layer_size]
			

	Return: RMSD_thresh_vals --- Dataframe containing the thresholds for bins and classes. 
	"""
	RMSD_thresh_vals={}
	for RMSD_mode, params in zip(mode, staggered_params):
		[staggering_value, max_angstrom, first_layer_size, second_layer_size]=params
		second_layer_list_down=np.arange(0,max_angstrom,second_layer_size) # makes 0,1,2,...4
		second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
		second_layer_list_up=np.append(second_layer_list_up,100)
		first_layer_thresh_list=[]

		for (up,down) in zip(second_layer_list_up,second_layer_list_down):
			staggerbins_down=np.arange(down,up,staggering_value)
			bin_list=[]
			x=first_layer_size
			for i in staggerbins_down:
				if x<up:
					x=i+first_layer_size
					x=np.around(x,2)
					i=np.around(i,2)
					bin_list.append((i,x))
			first_layer_thresh_list.append(bin_list)
		first_layer_thresh_list_down=np.arange(0,(max_angstrom-first_layer_size), staggering_value )
		first_layer_thresh_list_down=[np.around(x,2) for x in first_layer_thresh_list_down]
		first_layer_thresh_list_up=np.arange(first_layer_size, max_angstrom, staggering_value)
		first_layer_thresh_list_up=[np.around(x,2) for x in first_layer_thresh_list_up]
		first_layer_thresh_list=list(zip(first_layer_thresh_list_down,first_layer_thresh_list_up))

		RMSD_thresh_vals[RMSD_mode+"_threshold_list"]=first_layer_thresh_list_up
		RMSD_thresh_vals[RMSD_mode+"_threshup"]=second_layer_list_up
		RMSD_thresh_vals[RMSD_mode+"_threshlow"]=second_layer_list_down
		print(RMSD_thresh_vals)
	return RMSD_thresh_vals

def get_bin_nom_dic_equal(equal_params,mode):
	"""
	get equally spaced bins within min and max.
	step-size gives the number of bins combined to form a class. 
	"""
	RMSD_thresh_vals={}
	print(mode)
	print(equal_params)
	for RMSD_mode, equal_params_RMSD in zip(mode, equal_params):
		print(RMSD_mode)
		print(equal_params_RMSD)
		start_bins_staggered,stop_bins_staggered, nr_bins_staggered, bins_in_class_RMSD=equal_params_RMSD
		nr_bins_staggered=int(nr_bins_staggered)
		bins_in_class_RMSD=int(bins_in_class_RMSD)
		threshold_list=np.linspace(start_bins_staggered,stop_bins_staggered, nr_bins_staggered)
		
		threshold_list_nom=threshold_list[::bins_in_class_RMSD]
		print(threshold_list_nom)
		print(threshold_list)
		threshup=[]
		if threshold_list_nom[0]!=0:
			new_thresh=[0]
			threshlow=np.append(new_thresh, threshold_list_nom)
			threshup=threshlow[1]
		else:
			threshlow=threshold_list_nom


		#threshlow=np.append(threshlow,[100])
		threshup=np.append(threshup,threshold_list_nom[1:])
		threshup=np.append(threshup,[100])
		RMSD_thresh_vals[RMSD_mode+"_threshold_list"]=threshold_list
		RMSD_thresh_vals[RMSD_mode+"_threshup"]=threshup
		RMSD_thresh_vals[RMSD_mode+"_threshlow"]=threshlow
		print(RMSD_thresh_vals)
	return RMSD_thresh_vals

def make_RMSD_num_df(model_dir, actual_dir, mode):
	RMSD_num_df=pd.DataFrame()
	print(mode)
	for file in os.listdir(model_dir):
		file=os.path.join(model_dir, file)
		filename=os.path.splitext(os.path.basename(file))[0]
		filename=filename.replace(".pdb", "")
		actual_file=os.path.join(actual_dir, filename+".pdb")
		RMSD_num_row=get_RMSD_num(file, actual_file, mode)
		ID_df=pd.DataFrame()
		ID_df["ID"]=[filename]
		#RMSD_num_row = RMSD_num_row.fillna(value=np.nan)
		RMSD_num_row=pd.concat([RMSD_num_row, ID_df],axis=1)
		RMSD_num_df=pd.concat([RMSD_num_df, RMSD_num_row])
		RMSD_num_df.reset_index(drop=True,inplace=True)
	return RMSD_num_df

def make_val_df(balanced_params, staggered_params,equal_params,RMSD_num_df, mode):
	RMSD_num_df=None
	val_df=pd.DataFrame()
	if balanced_params is not None:
		print("balanced")
		RMSD_thresh_vals=get_bin_nom_dic_balanced(balanced_params,RMSD_num_df, mode)
	if staggered_params is not None:
		print("staggered")
		RMSD_thresh_vals=get_bin_nom_dic_staggered(staggered_params,mode)
	if equal_params is not None:
		print("equal")
		RMSD_thresh_vals=get_bin_nom_dic_equal(equal_params, mode)
	print(RMSD_thresh_vals)
	
	for RMSD_mode in mode:
		threshold_list=RMSD_thresh_vals[RMSD_mode+"_threshold_list"]
		threshup=RMSD_thresh_vals[RMSD_mode+"_threshup"]
		threshlow=RMSD_thresh_vals[RMSD_mode+"_threshlow"]

		threshold_list=np.asarray([round(x,2) for x in threshold_list])
		threshup=np.asarray([round(x,2) for x in threshup])
		threshlow=np.asarray([round(x,2) for x in threshlow])
		nom_values=np.arange(0,len(threshup))
		row = pd.Series({"nom_label":nom_values, "lower_threshold":threshlow, "upper_threshold":threshup, "bin_thresholds":threshold_list},name=RMSD_mode)
		val_df = val_df.append(row)
	
	return val_df


def get_full_dataset_df(model_dir, log_dir,actual_dir,mode, equal_params=None, balanced_params=None,staggered_params=None,save_to_path=None):
							
	"""
	This function is for transformaing all the individual feature and target Y dataframes of 
	each file into one big dataframe. Each row will represent one model file. 

	Input:  model_dir ---  Path of directory where all the model pdb files are saved that are to be
							used for training the model
			log_dir  ---   Path of directory containing all the abymod log files
			actual_dir --- Path of directory containing all 'real' pdb files taken from the AbDb database
			threshold_list ---  list of bin threshold values in Angstroms. Is the 
								RMSD is above or exactely this value, the bin value 
								will be 1, if the RMSd is 0. 
			threshlist --- list of tuples containing threshlow and threshup 
							(threshlow ---  List of lower threshold values of the nominal categories)
							(threshup ---  List of upper threshold values of the nominal categories)

	Return: full_dataset_df  --- Dataframe containing all model files with information of all features 
								that can be calculated using qualiloop library
							

	"""
	full_dataset_df=pd.DataFrame()
	RMSD_num_df=make_RMSD_num_df(model_dir, actual_dir, mode)

	#get bin and nominal cut-off values in the form of a df named val_df:
	val_df=make_val_df(balanced_params, staggered_params,equal_params,RMSD_num_df, mode)

	all_RMSD_df=pd.DataFrame()
	for RMSD_mode in mode:
		val_df_mode=val_df.loc[[RMSD_mode]]
		RMSD_num_df_mode=RMSD_num_df[[RMSD_mode,"ID"]]
		bin_thresholds=val_df_mode["bin_thresholds"]
		upper_threshold=val_df_mode["upper_threshold"]
		lower_threshold=val_df_mode["lower_threshold"]
		try:
			bin_thresholds=bin_thresholds.tolist()[0]
			upper_threshold=upper_threshold.tolist()[0]
			lower_threshold=lower_threshold.tolist()[0]
		except:
			pass
		nombin_col=make_df_nombin(RMSD_mode, None, None,bin_thresholds, upper_threshold, lower_threshold, RMSD_num_df=RMSD_num_df_mode)
		try:
			all_RMSD_df=all_RMSD_df.merge(nombin_col, on=['ID'])
		except:
			all_RMSD_df=pd.concat([nombin_col, all_RMSD_df])

	all_features_full=pd.DataFrame()
	for file in os.listdir(model_dir):
		file=os.path.join(model_dir, file)
		filename=os.path.splitext(os.path.basename(file))[0]
		filename=filename.replace(".pdb", "")
		actual_file=os.path.join(actual_dir, filename+".pdb")
		log_name=os.path.join(log_dir, filename+".log")
		if os.path.isfile(file)==False or os.path.isfile(actual_file)==False or os.path.isfile(log_name)==False:
			continue
		elif os.path.getsize(file)==0 or os.path.getsize(actual_file)==0 or os.path.getsize(log_name)==0:
			continue
		all_features_df=more_features(file,log_name)
		try:
			del all_features["index"]
		except:
			pass
		all_features_full=pd.concat([all_features_full, all_features_df])

			
	#all_features_df.reset_index(drop=True, inplace=True)
	#all_RMSD_df.reset_index(drop=True, inplace=True)
	#full_dataset_df.reset_index(drop=True, inplace=True)
	full_dataset_df=all_RMSD_df.merge(all_features_full, on=['ID'])
	full_dataset_df=full_dataset_df.merge(RMSD_num_df, on=['ID'])

	print("FINNNAAL")
	print(full_dataset_df)
	print(full_dataset_df.columns)


	#save df to csv
	if save_to_path==None:
		full_dataset_df.to_csv("full_dataset_df.csv", index=None)
	else:
		full_dataset_df.to_csv(save_to_path,index=None)

	return full_dataset_df, val_df


###############################################################################################################################################################
#Functions for building models and evaluating their performance

def split_full_dataset_df(full_dataset_df, target_y, del_y):
	"""
	Do a 1/3 to 2/3 split of the training data. This is done prior to Cross-validation,
	to create a separate validation dataset that was not used for triaining the model.
	y is a dataframe containing only the values of the target y (i.e. class to be predicted).
	X is a dataframe contianing all the features, but not any of the target y columns. X is 
	the input data for the model. 


	Input:  full_dataset_df  --- Dataframe containing all model files with information of all features 
								that can be calculated using qualiloop library
			target_y ---  The RMSD class to be predicted by the model
			del_y --- list of all column names which are RMSD class values 

	Return: X_train, X_test, y_train, y_test --- The train and test(more correctly validation here) dataframes of X and y
	"""
	y=full_dataset_df[target_y]
	X=full_dataset_df.drop(columns=del_y, errors='ignore')
	timeout_start = time.time()
	timeout=30 # max of 200 seconds
	condition=True
	while time.time() < timeout_start + timeout:
		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
		for a in [X_train, X_test, y_train, y_test]:
			if len(np.unique(a.to_numpy()))==1 :
				condition=False
		if len(np.unique(y_train))!=len(np.unique(y_test)): #make sure that all classes in teh train set are also represented in the test set
			condition=False		
		if condition:
			break
	if condition==False:
		return pd.DataFrame(),pd.DataFrame(),pd.DataFrame(),pd.DataFrame()

	return X_train, X_test, y_train, y_test

"""def define_RMSD_modes(full_dataset_df):
	
	Define lists of target ys corresponding to RMSD modes and list of all target ys. 


	Input:  full_dataset_df --- full dataset containing all features and target ys

	Return: del_y  ---  list of all column names that correspond to a target y
			local_CA_y_col_list  ---  list of all column names that correspond to a target y for local Calpha RMSD calculations
			local_AA_y_col_list  ---  list of all column names that correspond to a target y for local atom RMSD calculations
			global_CA_y_col_list  ---  list of all column names that correspond to a target y for global Calpha RMSD calculations
			global_AA_y_col_list  ---  list of all column names that correspond to a target y for global atom RMSD calculations
			numy_y_col  ---  list of the numercial target y columns
	
	#make lists of Y-values for local_CA
	local_CA_bin_cols = [col for col in full_dataset_df if col.startswith('local_CA_bin')]
	local_CA_nom_cols= [col for col in full_dataset_df if col.startswith('local_CA_nom')]
	local_CA_y_col_list=local_CA_bin_cols+local_CA_nom_cols

	#make lists of Y-values for local_AA
	local_AA_bin_cols = [col for col in full_dataset_df if col.startswith('local_AA_bin')]
	local_AA_nom_cols= [col for col in full_dataset_df if col.startswith('local_AA_nom')]
	local_AA_y_col_list=local_AA_bin_cols+local_AA_nom_cols
	#make lists of Y-values for global_CA
	global_CA_bin_cols = [col for col in full_dataset_df if col.startswith('global_CA_bin')]
	global_CA_nom_cols= [col for col in full_dataset_df if col.startswith('global_CA_nom')]
	global_CA_y_col_list=global_CA_bin_cols+global_CA_nom_cols
	#make lists of Y-values for global_AA
	global_AA_bin_cols = [col for col in full_dataset_df if col.startswith('global_AA_bin')]
	global_AA_nom_cols= [col for col in full_dataset_df if col.startswith('global_AA_nom')]
	global_AA_y_col_list=global_AA_bin_cols+global_AA_nom_cols

	#make list of numeric y columns
	num_y_col=["local_CA","local_AA", "global_CA", "global_AA"]

	#delete all other y-columns from the X_df by making a list of all other Y-values and droping them
	del_y=local_CA_y_col_list+ local_AA_y_col_list+global_CA_y_col_list+ global_AA_y_col_list +num_y_col

	return del_y, local_CA_y_col_list, local_AA_y_col_list,global_CA_y_col_list, global_AA_y_col_list, num_y_col

"""
def clean_dataset(df):
	"""
	Clean the input dataframe. Checks it is dataframe, removes any empty rows (i.e any model structures for which a 
	feature calculation gave no result) or rows with nan, infinity or -inifity values. 
	Returns the cleaned dataframe. 


	Input:  df  --- Dataframe to be cleaned

	Return: cleaned dataframe
	"""
	assert isinstance(df, pd.DataFrame)
	#df.dropna(inplace=True)
	indices_to_keep = ~df.isin([np.nan, np.inf, -np.inf, "nan", "None"]).all(1)
	return df

def all_vals_same(column):
	"""
	Check if the input column only has the same values. This check is implimented so that a target y 
	doesnt only have instances of one class in the dataset, which would falsify the performance measurements


	Input:  column  --- Dataframe column to be checked

	Return: True or False (True is all values are the same, otherwise False)
	"""
	if len(np.unique(column)) == 1:
		return True
	else: 
		return False



def all_classifier_types(model_list=None):
	"""
	Define all classifiers that will be tested. This set can be changed using the model_list input.  
	IMPORTANT: For Hyperparameter optimization the name of the model has to be set to the official 
	scikit name, otherwise the program hyperparam_opti.py will not recognize the model type. 

	If a model_list is added, all models must also be found in hyperparam_opti.py!!!


	Input:  model_list (optional) --- list of model names and models (same format as seen below), if
										not default of none, only the models in this list are used to 
										build the model
		
	Return: models --- list of model names and models to be evaluated and used for building the final model
	"""
	if model_list:
		return model_list
	models=[]
	models.append(('LogisticRegression', LogisticRegression(solver='liblinear', multi_class='ovr')))
	#models.append(('LinearDiscriminantAnalysis',LinearDiscriminantAnalysis()))
	models.append(('KNeighborsClassifier', KNeighborsClassifier()))
	models.append(('DecisionTreeClassifier', DecisionTreeClassifier()))
	models.append(('GaussianNB', GaussianNB()))
	models.append(('RandomForestClassifier',RandomForestClassifier(n_estimators=200)))
	no_vote_models=models.copy()
	models.append(('SVC', SVC(gamma='auto')))
	no_vote_models.append(('SVC', SVC(gamma='auto', probability=True)))  #without defining probability=True, soft voting fails
	models.append(('Voting_all_soft', VotingClassifier(estimators=no_vote_models, voting='soft', weights=[1]*len(no_vote_models), flatten_transform=True)))
	models.append(('XGBoost', xgb.XGBClassifier()))
	#models.append(('Voting_all_hard', VotingClassifier(estimators=no_vote_models, voting='hard'))) #has no proba prediction
	return models


def make_models(X_train, y_train, X_test, y_test,models,model_list=None, only_score=False):
	"""
	Build and train 7 different common classification models using the input X and y training set. 
	Each model is assessed on performance using K-fold Cross-Validation (10 folds). The average 
	Mathwes Correlation Coefficient of all folds for each model is stored in a dictionary. 


	Input:  X_train ---  Dataframe containing all features of the full dataset the model uses as input 
			y_train ---  Dataframe containing the target y values for each model structure in the dataset

	Return: MCC_dic ---  Dictionary containing two key-value pairs per model structure. 
							1) "Model_name":model name
							2) "MCC_mean": mean of all MCC values collected for 10 folds
	"""
	full_y=np.concatenate((np.unique(y_train), np.unique(y_test)))
	num_class=len(np.unique(full_y))
	if num_class<=1:
		num_class=2
	MCC_dic={"Model_name":[], "MCC_mean":[]}
	
	for (name, model) in models:
		print(model)
		MCC_per_fold=[]
		fold_no = 1
		strat_k_fold = RepeatedStratifiedKFold(n_splits=10,random_state=42,n_repeats=4)
		#k_fold = KFold(n_splits=10, shuffle=True)
		fold_dic={}
		#all_mccs = cross_val_score(estimator=model, X=X_train, y=y_train, cv=10, scoring=mcc_scorer)

		X_train.reset_index(drop=True, inplace=True)

		for Ktrain_index, Ktest_index in strat_k_fold.split(X_train, np.zeros(shape=(X_train.shape[0], 1))):
			KX_train, KX_test = X_train.iloc[Ktrain_index], X_train.iloc[Ktest_index]
			Ky_train, Ky_test = y_train.iloc[Ktrain_index], y_train.iloc[Ktest_index]
			if all_vals_same(Ky_train.to_numpy())==True:
				continue
			if name=="XGBoost" or name=="XGBClassifier":
				#KX_train, KX_test, Ky_train, Ky_test = KX_train.to_numpy(), KX_test.to_numpy(), Ky_train.to_numpy(), Ky_test.to_numpy()
				param={'objective':'binary:logistic'}
				if num_class>2:
					param['objective']='multi:softmax'
					param['num_class']=num_class+1
				model=model.set_params(**param)
				dMatrixTrain = xgb.DMatrix(KX_train, label=Ky_train)
				dMatrixTest = xgb.DMatrix(KX_test, label=Ky_test)
				if only_score==False:
					xgbT = xgb.train(param, dMatrixTrain, 10)
					rf_prediction = xgbT.predict(dMatrixTest)
				else:
					rf_prediction=model.predict(KX_test.to_numpy())
				print(rf_prediction)
				rf_prediction=rf_prediction.astype(np.float)
			else:
				if only_score==False:
					KX_train.reset_index(drop=True, inplace=True)
					Ky_train.reset_index(drop=True, inplace=True)
					model.fit(KX_train, Ky_train)
				try:
					rf_prediction = model.predict(KX_test)# Evaluations
				except:
					rf_prediction = model.predict(KX_test.values.flatten())# Evaluation
			try:
				MCC_of_fold=matthews_corrcoef(Ky_test,rf_prediction, sample_weight=None )
			except:
				rf_prediction = rf_prediction>0.5  # transorms into true, false -> as xgboost predicts a probability by default
				MCC_of_fold=matthews_corrcoef(Ky_test,rf_prediction, sample_weight=None )
			MCC_per_fold.append(MCC_of_fold)
			fold_dic[MCC_of_fold]=[fold_no]
			fold_no = fold_no + 1
		MCC_dic["Model_name"].append(str(name))
		MCC_dic["MCC_mean"].append(np.mean(MCC_per_fold))
	return MCC_dic

def hyperparameter_optimization(models, X_train, y_train, X_test, y_test):
	#uses hyperparam_opti.py
	#incorporate genetic algorithm here and other methods
	
	#models_opti=opti_gridsearch(models,X_train.to_numpy(), y_train.to_numpy(), X_test.to_numpy(), y_test.to_numpy())
	for i in range(5):
		try:
			models_opti=opti_genetic(models,X_train, y_train, X_test, y_test)
			MCC_dic_opti=make_models(X_train, y_train, X_test, y_test,models_opti)
			return models_opti, MCC_dic_opti
		except:
			pass
	MCC_dic=make_models(X_train, y_train, X_test, y_test,models)
	return models, MCC_dic

	
	


def get_best_submodel_opti(X_train, X_test, y_train, y_test,models, no_opti=False):
	"""
	Determine which of the models generated for a specfific target_y is the best 
	(which classifier, and if optimized or not). Returns a single-row dataframe 
	with information on model-type, model performance and what the target class is.
	Also contains the model itself. 
	Will later be assembled into a larger dataframe


	Input:  full_dataset_df ---  Dataframe containing all features of the full dataset and all target y's
			target_y ---  target class name to be used in this run
			del_y ---  list of the names of all other target y's

	Return: MCC_dic_opti ---  Dictionary containing two key-value pairs per model structure. 
							1) "Model_name":model name
							2) "MCC_mean": mean of all MCC values collected for 10 folds
			MCC_dic  --- same as MCC_dic_opti but with the non-optimized models
			models_torun  --- single-row dataframe with information on model-type, 
								model performance, what the target class is and the model itself.
	"""
	
	
	#as if X_train is empty, the split_full_dataset_df function could not create a train/test set that 
	#fulfills both the requirement, that not all the class labels are the same, and also that y-train and y-test 
	#contain the same number of class labels
	MCC_dic=make_models(X_train, y_train, X_test, y_test, models)
	print("MCC_DIC", MCC_dic)
	if no_opti==True:
		models_opti=models
		MCC_dic_opti=MCC_dic
	else:
		models_opti, MCC_dic_opti=hyperparameter_optimization(models, X_train, y_train, X_test, y_test)
	gc.collect()


	#--------
	#create best_models and MCC dataframes for analysis purposes. 
	if len(MCC_dic_opti.values()) >=1:
		MCC_df_opti=pd.DataFrame.from_dict(MCC_dic_opti, orient = 'columns')
		max_value_opti = max(MCC_df_opti["MCC_mean"].to_list())  # maximum value
		max_keys_opti=list(MCC_df_opti.loc[(MCC_df_opti['MCC_mean'] == max_value_opti, "Model_name")])

	if len(MCC_dic.values()) >=1:
		MCC_df=pd.DataFrame.from_dict(MCC_dic, orient = 'columns')
		max_value = max(MCC_df["MCC_mean"].to_list())  # maximum value
		max_keys=list(MCC_df.loc[(MCC_df['MCC_mean'] == max_value, "Model_name")])

	#------
	#select best model (if optimized or not) and add this to the models_torun dataframe.

	#visualizer.opti(MCC_df, MCC_df_non_optimized, models_torun)

	models_torun=pd.DataFrame()
	print(max_keys, max_keys_opti)
	if max_value_opti > max_value:
		for model_name, model_opti in models_opti:
			print(model_name, max_keys_opti)
			if [model_name]==max_keys_opti:
				models_torun["model_name"]= max_keys_opti
				models_torun["MCC"]= [max_value_opti]
		print("hhhhhhhopti",max_value_opti)
		print(X_train,y_train)
		#input(models_torun)
		try:
			X_train.drop(['index'], axis=1,inplace=True)
			X_test.drop(["index"],axis=1,inplace=True)
			try:
				X_train.drop(['index'], axis=1,inplace=True)
				X_test.drop(["index"],axis=1,inplace=True)
			except:
				pass
		except:
			pass
		model_opti.fit(X_train, y_train)
		return models_torun, MCC_df, MCC_df_opti, model_opti
	elif max_value_opti <= max_value:
		for model_name, model in models:
			if [model_name]==max_keys:
				models_torun["model_name"]= max_keys
				models_torun["MCC"]= [max_value]
		print("hhhhhhh",max_value)
		#input(models_torun)
		print(X_train,y_train)

		try:
			X_train.drop(['index'], axis=1,inplace=True)
			X_test.drop(["index"],axis=1,inplace=True)
		except:
			pass
		try:
			X_train.drop(['index'], axis=1,inplace=True)
			X_test.drop(["index"],axis=1,inplace=True)
		except:
			pass
		print("HERE")
		print(type(X_train),type(y_train))
		model.fit(X_train, y_train)
		return models_torun, MCC_df, MCC_df_opti, model



def first_layer(RMSD_mode, models,X_train, X_test, y_train_all, y_test_all, model_type, only_first=False, run_dir=None, no_opti=False):
	"""
	Create models for all binary bins. Save the predictions (in probability format) as a df to be used as 
	features in the second layer of the model. 


	Input:  full_dataset_df ---  Dataframe containing all features of the full dataset and all target y's

	Return: full_dataset_firstlayer_df  ---  Dataframe containing all features of the full dataset and all predictions for the bins
	"""
	
	MCC_of_bins=pd.DataFrame()
	prediction_df_train=pd.DataFrame()
	prediction_df_test=pd.DataFrame()
	if not os.path.exists(os.path.join(run_dir,model_type[0],"layer_1")):
			os.makedirs(os.path.join(run_dir,model_type[0],"layer_1"))
	for target_y in RMSD_mode:
		y_train=y_train_all[target_y]
		y_test=y_test_all[target_y]
		if X_train.empty or X_test.empty or y_train.empty or y_test.empty:
			continue
		if "nom" in target_y:
			continue
		if os.path.exists(os.path.join(run_dir,model_type[0],"layer_1","bin"+str(target_y)+".pkl")):
			best_model_of_bin=joblib.load(os.path.join(run_dir,model_type[0],"layer_1","bin"+str(target_y)+".pkl"))
		elif os.path.exists(os.path.join(run_dir,model_type[0],"layer_1","bin"+str(target_y)+".json")):
			best_model_of_bin=xgb.XGBClassifier()
			best_model_of_bin.load_model(os.path.join(run_dir,model_type[0],"layer_1","bin"+str(target_y)+".json"))
		else:
			print(y_test)
			if no_opti==True:
				models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models, no_opti=True)
			else:
				models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models)
			if not os.path.isdir(os.path.join(run_dir,model_type[0],"layer_1")):
				os.makedirs(os.path.join(run_dir,model_type[0],"layer_1"))

			if type(best_model_of_bin).__name__=="XGBClassifier":
				best_model_of_bin.save_model(os.path.join(run_dir,model_type[0],"layer_1","bin"+str(target_y)+".json"))
			else:
				joblib.dump(best_model_of_bin, os.path.join(run_dir,model_type[0],"layer_1","bin"+str(target_y)+".pkl")) 

		prediction_train = best_model_of_bin.predict_proba(X_train)
		print("First_LAYER_INPUT",X_train.columns)
		prediction_df_train[target_y]=prediction_train[:,1]
		#MCC_target=make_models(X_train, y_train, X_test, y_test,[(type(best_model_of_bin).__name__, best_model_of_bin)],only_score=True)
		#MCC_val=MCC_target["MCC_mean"]
		#MCC_of_bins[target_y]=MCC_val

		prediction_test= best_model_of_bin.predict_proba(X_test)
		prediction_df_test[target_y]=prediction_test[:,1]

	return prediction_df_train, prediction_df_test, MCC_of_bins

def staggered_second_layer(X_train, X_test, y_train_all, y_test_all, prediction_df_train,prediction_df_test,MCC_of_bins,val_df, RMSD_mode, models,model_type, factor_MCC=False, run_dir=None,no_opti=False):
	"""
	Create input dataframe containing all features and transformed values for the predicted bin probabilities calculated in 
	the first layer. The probabilities for bin1, bin1.1,... are multiplied by the MCC of the model that predicted the value and these
	products are then summed for each nominal category. This sum is divided by the number of bins within the nominal category.
	I.e for nom1, which is 1-2, The value will be:
	(bin1*MCC1+bin1.1*MCC1.1+,...+bin1.9*MCC )/10

	Input:  full_dataset_df ---  Dataframe containing all features of the full dataset and all target y's

	Return: full_dataset_firstlayer_df  ---  Dataframe containing all features of the full dataset and all target y's
	"""
	gc.collect()
	print(val_df)
	threshlow=val_df[["lower_threshold"]]
	threshup=val_df[["upper_threshold"]]
	for target_y in RMSD_mode:
		if "nom" in target_y:
			break
   
	print(prediction_df_train)
	print(prediction_df_test)
	for train_test in ("prediction_df_train","prediction_df_test"):
		second_layer_df=pd.DataFrame()
		print("STRT:",train_test)
		if train_test=="prediction_df_train":
			in_use=prediction_df_train
		else:
			in_use=prediction_df_test
		for i, i_low in zip(threshup[0], threshlow[0]):
			print("THRESHUP;THRESHLOW",i,i_low)
			if i==100:
				continue
			else:
				multiply_list=[]
				for key,value in in_use.items():
					key=key.split("bin", 1)[1]
					key=float(key)
					if key <=i and key>i_low:
						#if factor_MCC:
						#	value=value*MCC_of_bins[key].astype(np.float)
						multiply_list.append([value])
				number_of_vals=len(multiply_list)
				arr = np.array(multiply_list)
				sum_array=arr.sum(axis=0)/number_of_vals
				sum_array=sum_array[0]   
				print("SUMM ARRAY", sum_array)
				print(second_layer_df)
				print("ISSUE HERE",second_layer_df.shape, sum_array.shape)
				second_layer_df[str("nom")+str(i)]=sum_array
		if train_test=="prediction_df_train":
			second_layer_df_train=second_layer_df
		else:
			second_layer_df_test=second_layer_df  
	print("second_layer_df_train",second_layer_df_train)
	print("second_layer_df_test",second_layer_df_test)
	X_train.reset_index(drop=True, inplace=True)
	second_layer_df_train.reset_index(drop=True, inplace=True)
	second_layer_df_test.reset_index(drop=True, inplace=True)
	X_test.reset_index(drop=True, inplace=True)
	X_train_with_preds=pd.concat([X_train, second_layer_df_train],axis=1)
	X_train_with_preds.to_csv(os.path.join(run_dir,"second_layer_df_train.txt"), index=None)

	X_test_with_preds=pd.concat([X_test, second_layer_df_test],axis=1)
	X_test_with_preds.to_csv(os.path.join(run_dir,"second_layer_df_test.txt"), index=None)
	X_test_with_preds.reset_index(drop=True, inplace=True)
	X_train_with_preds.reset_index(drop=True, inplace=True)
	print(X_train_with_preds.shape, X_test_with_preds.shape)

	
	y_train=y_train_all[target_y]
	y_test=y_test_all[target_y]
	X_train=X_train_with_preds
	X_test=X_test_with_preds

	print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
	print(X_train.columns, X_test.columns)
	X_train.reset_index(drop=True, inplace=True)
	X_test.reset_index(drop=True, inplace=True)
	y_train.reset_index(drop=True, inplace=True)
	y_test.reset_index(drop=True, inplace=True)
	print(y_test)
	print(y_train)
	if not os.path.exists(os.path.join(run_dir)):
			os.makedirs(os.path.join(run_dir))
	try:
		X_train.drop(['index'], axis=1,inplace=True)
		X_test.drop(["index"],axis=1,inplace=True)
		try:
			X_train.drop(['index'], axis=1,inplace=True)
			X_test.drop(["index"],axis=1,inplace=True)
		except:
			pass
	except:
		pass

	if no_opti==True:
		models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models, no_opti=True)
	else:
		models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models)

	if type(best_model_of_bin).__name__=="XGBClassifier":
		best_model_of_bin.save_model(os.path.join(run_dir,model_type[0],"all"+".json"))
	else:
		joblib.dump(best_model_of_bin, os.path.join(run_dir,model_type[0],"all"+".pkl")) 
	best_model_MCC_secondlayer=models_torun["MCC"]
	best_model_of_bin_all=best_model_of_bin
	models_torun_all=models_torun
	print("second_layer MCC:")
	print(best_model_MCC_secondlayer)
		

	y_train=y_train_all[target_y]
	y_test=y_test_all[target_y]
	X_train=second_layer_df_train
	X_test=second_layer_df_test

	print(X_train.shape, X_test.shape, y_train.shape, y_test.shape)
	print(X_train.columns, X_test.columns)
	X_train.reset_index(drop=True, inplace=True)
	X_test.reset_index(drop=True, inplace=True)
	y_train.reset_index(drop=True, inplace=True)
	y_test.reset_index(drop=True, inplace=True)
	print(y_test)
	print(y_train)
	try:
		X_train.drop(['index'], axis=1,inplace=True)
		X_test.drop(["index"],axis=1,inplace=True)
		try:
			X_train.drop(['index'], axis=1,inplace=True)
			X_test.drop(["index"],axis=1,inplace=True)
		except:
			pass
	except:
		pass
	if no_opti==True:
		models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models, no_opti=True)
	else:
		models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models)
	if type(best_model_of_bin).__name__=="XGBClassifier":
		best_model_of_bin.save_model(os.path.join(run_dir,model_type[0],"only_sums"+".json"))
	else:
		joblib.dump(best_model_of_bin, os.path.join(run_dir,model_type[0],"only_sums"+".pkl")) 
	best_model_MCC_secondlayer=models_torun["MCC"]
	best_model_of_bin_only_sums=best_model_of_bin
	models_torun_only_sums=models_torun
	
	print("second_layer MCC:")
	print(best_model_MCC_secondlayer)

	return models_torun_all, best_model_of_bin_all, models_torun_only_sums,best_model_of_bin_only_sums,X_train_with_preds,X_test_with_preds,second_layer_df_train,second_layer_df_test

	

def double_staggered_run(X,y, del_y,models, RMSD_mode, val_df, model_type,run_dir=None, no_opti=False):
	"""
	Full run for the double staggered model (see documentation for explanation).  

	Input:  full_dataset_df ---  Dataframe containing all features of the full dataset and all target y's
	
	"""
	#for visualizer
	MCC_all=pd.DataFrame()
	MCC_all_non_optimized=pd.DataFrame()

	print(RMSD_mode)
	for target_y in RMSD_mode:
		if "nom" in target_y:
			X_train, X_test, y_train_all, y_test_all = train_test_split(X, y, stratify=y[target_y], test_size=0.2)
			print(X_train.shape, X_test.shape, y_train_all.shape, y_test_all.shape)
			

	prediction_df_train, prediction_df_test, MCC_of_bins=first_layer(RMSD_mode, models,X_train, X_test, y_train_all, y_test_all,
															model_type, run_dir=run_dir, no_opti=no_opti)

	(models_torun, model, models_torun_only_sums, model_only_sums,X_train_with_preds,X_test_with_preds,second_layer_df_train,second_layer_df_test)=staggered_second_layer(X_train, X_test, y_train_all, y_test_all, 
																					prediction_df_train,prediction_df_test,MCC_of_bins,val_df, RMSD_mode, 
																					models, model_type, run_dir=run_dir, no_opti=no_opti)
	y_train=y_train_all[target_y]
	y_test=y_test_all[target_y]
	return y_train,y_test,models_torun, model, models_torun_only_sums, model_only_sums,MCC_of_bins,X_train_with_preds,X_test_with_preds,second_layer_df_train,second_layer_df_test


if __name__ == '__main__':
	print("START")
	#full_experimentation()
	model_dir=sys.argv[1]
	log_dir=sys.argv[2]
	actual_dir=sys.argv[3]
	full_RMSD_df=pd.DataFrame()
	

	for file in os.listdir(model_dir):
		filename=file.replace(".pdb.model", "")
		file=os.path.join(model_dir, file)
		actual_file=os.path.join(actual_dir, filename+".pdb")
		print(filename)
		print(actual_file)
		RMSD_num_df=get_RMSD_num(file, actual_file,mode)
		ID_df=pd.DataFrame()
		ID_df["ID"]=[filename]
		print(ID_df)
		RMSD_num_df = pd.concat([RMSD_num_df, ID_df], axis=1)
		print(RMSD_num_df)
		full_RMSD_df = pd.concat([RMSD_num_df, full_RMSD_df], axis=0)
	full_RMSD_df.to_csv("~/sync_project/AMA_I/RMSD_vals_"+"AMA_I.csv")


	"""import seaborn as sns
	global_CA = full_RMSD_df["global_CA"].to_list()
	global_CA=[float(x) for x in global_CA]
	full_RMSD_df["global_CA"]=global_CA
	ax = sns.barplot(x="ID", y="global_CA", data=full_RMSD_df)
	ax.set_title('RMSD values of the CDR-H3 region for abYmod - AMA-I')
	ax.set(xlabel="PDB code", ylabel = 'RMSD (Å)')
	ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
	fig=ax.get_figure()
	fig.savefig(os.path.join("/home/lilian/sync_project/AMA_I/","bar_chart_AMA_I.png"))"""

	"""model_dir=sys.argv[1]
	log_dir=sys.argv[2]
	actual_dir=sys.argv[3]
	equal_params=[[1,8,8], [1,8,8], [1,8,8],[1,8,8]]
	bins_in_class=[2,2,2,2]
	nr_of_balanced_bins=[10,10,10,10]
	#full_dataset_df, val_df=get_full_dataset_df(model_dir, log_dir,actual_dir,equal_params=i, 
	#	bins_in_class=bins_in_class,save_to_path=os.path.join("./new_full_df"+str(i)),
	#	nr_of_balanced_bins=nr_of_balanced_bins, balanced=False)
	full_dataset_df, val_df=get_full_dataset_df(model_dir, log_dir,actual_dir,equal_params=equal_params, 
	 bins_in_class=bins_in_class,save_to_path="./full_dataset_balanced10bins_10_04.csv",
	 nr_of_balanced_bins=nr_of_balanced_bins, balanced=True)

	"""
	"""full_dataset_df=pd.read_csv(sys.argv[4],header=0)
	best_models, MCC_df, best_models_non_optimized, MCC_df_non_optimized=full_run(full_dataset_df)
	best_models.to_csv("./best_models.csv", index=None)
	MCC_df.to_csv("./MCCs.csv", index=None)
	best_models_non_optimized.to_csv("./best_models_non_opti.csv", index=None)
	MCC_df_non_optimized.to_csv("./MCCs_non_opti.csv", index=None)
	print(best_models)
	print(MCC_df)"""