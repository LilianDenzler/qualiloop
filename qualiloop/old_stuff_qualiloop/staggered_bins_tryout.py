#!/usr/bin/env python
import sys
import os
import pandas as pd
import math
from pathlib import Path
import numpy as np

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import style

# Evaluations
from sklearn.metrics import classification_report,confusion_matrix
# Random Forest
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
from tensorflow.keras.models import Sequential, save_model
import joblib
from sklearn.preprocessing import MinMaxScaler
import collections
import shutil
from shutil import copyfile
'''original setting:
	staggering_value=0.1
	max_angstrom=5
	first_layer_size=0.5
	second_layer_size=1
	'''

def RMSD_binary(threshold,RMSD_df):
	get_all_bins(RMSD_mode, file,actual_file, threshold_list, RMSD_num_df=None, fast_mode=None)

	return (dfnew)

def RMSD_nom(RMSD_df,second_layer_list_up,second_layer_list_down):
	modelmaker.RMSD_nom(RMSD_mode, threshlow,threshup, file, actual_file,counter=None, RMSD_num_df=None, fast_mode=None)
	return (dfnew)


def model(X,y, bin_models,up):
	get_best_submodel_opti(X_train, X_test, y_train, y_test,models)

	return np.mean(MCC_per_fold)

def make_bin_df(full_dataset_df,bin_models,staggering_value,max_angstrom,first_layer_size,second_layer_size,accuracy_csv_path, features_to_include):
	#here the staggering is set
	#total range: 0-5 i.e above 5Angstrom there is no more differentiation
	#first_layer bins are 1 Angstrom in size, second layer bins are 2Angstroms in size
	#staggering value determines the shift from bin to bin in the first layer
	second_layer_list_down=np.arange(0,max_angstrom,second_layer_size) # makes 0,1,2,...4
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
	second_layer_list_up=np.append(second_layer_list_up,max_angstrom)
	first_layer_thresh_list=[]

	df=RMSD_nom(full_dataset_df,list(second_layer_list_up),list(second_layer_list_down)) # gives df of nominal values
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
	first_layer_thresh_list=zip(first_layer_thresh_list_down,first_layer_thresh_list_up)

	for down,up in first_layer_thresh_list:
		df2=RMSD_binary(up,full_dataset_df)   #get binary value df

def make_bin_models(RMSD_name, first_layer_thresh_list):
	accuracy_bins={}
	for down,up in first_layer_thresh_list:
		print("_____________________ bins",down,up)

		models_torun, MCC_df, MCC_df_opti, best_model_of_bin=get_best_submodel_opti(X_train, X_test, y_train, y_test,models)
		if not os.path.isdir(os.path.join(run_dir,"bin_models","layer_1")):
			os.makedirs(os.path.join(run_dir,"bin_models","layer_1"))
		joblib.dump(best_model_of_bin, os.path.join(run_dir,"bin_models","layer_1","bin"+str(target_y)+".pkl"))

	return df  #df of nominal values

def make_test_df(bin_models, test_file,second_layer_size, max_angstrom,first_layer_size,staggering_value, features_to_include): #accuracy_csv_path
	#####test_df2
	predicted_vals={}
	for file in os.listdir(bin_models):
		filename=file.replace(".pkl","")
		filename=filename.replace("bin","")
		model= joblib.load(os.path.join(bin_models,file))
		vals=test_df2.iloc[:,:].to_numpy()
		predictions=model.predict(vals)
		predicted_vals.update({filename:predictions})
	
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
	second_layer_list_up=np.append(second_layer_list_up,5)
	summed_layer_two={}
	#accuracy_dic=MCC_weighting(accuracy_bins,second_layer_list_up, predicted_vals, first_layer_size,staggering_value,max_angstrom)
	for i in second_layer_list_up:
		i_low=i-1
		multiply_list=[]
		accuracies=[]
		for key,value in predicted_vals.items():
			key=float(key)
			if key <=i and key>i_low:
				#multiply_by=accuracy_dic[key]
				multiply_list+=[value] #*multiply_by
				#print("HEEEEEEEEEEEEEY",type(value))
		number_of_vals=len(multiply_list)
		arr = np.array(multiply_list)
		sum_array=arr.sum(axis=0)/number_of_vals   
		#print("sum array",sum_array)
		#sum_array = sum_array/number_of_vals
		summed_layer_two.update({i:sum_array})
	test_df=pd.DataFrame.from_dict(summed_layer_two)
	print("test_df",test_df)
	test_df.columns = test_df.columns.astype(str)
	return test_df



def third_layer(a, third_model_path):
	a.isnull().any()
	a_nom=a[["local_CA_nom"]].to_numpy()
	a_vals=a.drop(labels="local_CA_nom",axis=1).to_numpy()
	up="THIRD_LAYER"
	model(a_vals, a_nom, third_model_path, up)



def run(RMSD_file, feature_csv, bin_models, test_file,feature_directory,accuracy_csv_path,third_model_path):
	staggering_value=0.2
	max_angstrom=5
	first_layer_size=1
	second_layer_size=2
	
	#make and dump binary models
	nom_df=make_bin_df(RMSD_file, feature_csv, bin_models,staggering_value,max_angstrom,first_layer_size,second_layer_size, accuracy_csv_path, features_to_include_noID)
	#make prediction_df
	test_df=make_test_df(bin_models, test_file, second_layer_size,max_angstrom,first_layer_size,staggering_value,  features_to_include_noID)
	#join prediction df to feature csv
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size)
	second_layer_list_up=[str(x) for x in second_layer_list_up]
	second_layer_list_up.append(str(max_angstrom))
	
	list_features=second_layer_list_up+features_to_include
	a=a[list_features]
	a=a.merge(nom_df,on="ID")
	a=a.drop(labels=["ID","local_AA_nom","global_AA_nom","global_CA_nom"], axis=1)
	
	third_layer(a, third_model_path)
	a.to_csv(os.path.join(feature_directory,"MERGED_STAGGERED_tmp"), index=False)








run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[5])
#run_with_unknown(sys.argv[1],sys.argv[2])
#run_with_unknown(sys.argv[3],sys.argv[4])



#with setting:

'''staggering_value=0.5
	max_angstrom=5
	first_layer_size=1
	second_layer_size=2'''
#####=========0.6266835110483702

'''staggering_value=0.2
	max_angstrom=5
	first_layer_size=1
	second_layer_size=2'''
#####=========MCC=0.6471775637524335 		!!!!

'''staggering_value=0.2
	max_angstrom=5
	first_layer_size=1.5
	second_layer_size=3'''
####=========0.6274200759280766
######################################################################################

'''staggering_value=0.2
	max_angstrom=5
	first_layer_size=1.5
	second_layer_size=2'''
####MCC=0.6572073279745841





#python3 staggered_bins.py ~/sync_project/Feature/train_staggered.csv ~/sync_project/Feature/train_staggered.csv ~/sync_project/staggered_models ~/sync_project/Feature/train_staggered.csv ~/sync_project/Feature/ ~/sync_project/Feature/accuracy_df.csv ~/sync_project/Feature/test_staggered.csv ~/sync_project/Feature/binTHIRD_LAYER.pkl