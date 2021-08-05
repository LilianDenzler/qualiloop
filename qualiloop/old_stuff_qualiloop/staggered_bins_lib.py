#!/usr/bin/env python
import sys
import os
import pandas as pd
import math
import numpy as np

import pandas as pd


def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
#use this to bypass sci-kits forced warning when unpickling the model
import joblib

'''original setting:
	staggering_value=0.1
	max_angstrom=5
	first_layer_size=0.5
	second_layer_size=1
	'''


def RMSD_nom(RMSD_df,second_layer_list_up,second_layer_list_down):
	columns=["local_AA_nom","local_CA_nom","global_AA_nom", "global_CA_nom", "ID"]
	data=[]
	local_AA= RMSD_df.local_AA.tolist()
	local_CA= RMSD_df.local_CA.tolist()
	global_AA= RMSD_df.global_AA.tolist()
	global_CA=RMSD_df.global_CA.tolist()
	name=RMSD_df.ID.tolist()
	for (a,b,c,d,n) in zip(local_AA, local_CA, global_AA, global_CA, name):
		for (threshold_up,threshold_low) in zip(second_layer_list_up, second_layer_list_down):
			#print(threshold_low,threshold_up)
			#print(type(a))
			if threshold_up==max(second_layer_list_up):
				threshold_up=50
			if type(a)==str:
				pass
			elif a>=threshold_low and a<threshold_up:
				a=">"+str(threshold_low)+"<"+str(threshold_up)
			if type(b)==str:
				pass
			elif b>=threshold_low and b<threshold_up:
				b=">"+str(threshold_low)+"<"+str(threshold_up)
			if type(c)==str:
				pass
			elif c>=threshold_low and c<threshold_up:
				c=">"+str(threshold_low)+"<"+str(threshold_up)
			if type(d)==str:
				pass
			elif d>=threshold_low and d<threshold_up:
				d=">"+str(threshold_low)+"<"+str(threshold_up)
		write=[a,b,c,d,n]
		#print(write)
		data.append(write)
	dfnew = pd.DataFrame(data, columns=columns)
	#print(dfnew)
	return (dfnew)


def make_test_df(bin_models, test_df,second_layer_size, max_angstrom,first_layer_size,staggering_value, features_to_include): #accuracy_csv_path
	test_df.isnull().any()
	#accuracy_bins=pd.read_csv(accuracy_csv_path, header=0)
	#accuracy_bins=accuracy_bins.to_dict()

	test_df2 = test_df.copy()
	test_df2 = pd.DataFrame(data=test_df2, columns=features_to_include)
	#print(test_df2.head())
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
	#print("test_df",test_df)
	test_df.columns = test_df.columns.astype(str)
	return test_df


def run_with_unknown(unknown_df, third_model_path,bin_model_path,test=None) :
	staggering_value=0.2
	max_angstrom=5
	first_layer_size=1
	second_layer_size=2
	features_to_include=['tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff'] #deleted tip_res
	second_layer_list_down=np.arange(0,max_angstrom,second_layer_size) # makes 0,1,2,...4
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size) # makes 1,2,...4
	second_layer_list_up=np.append(second_layer_list_up,max_angstrom)
	bin_models=os.path.join(bin_model_path)
	#first we take our csv with all our features and remove all unwanted features (and calculate local_CA_nom class if you want for testing)

	#we then feed into first layer and calculate the prediction values for each categories from second layer back into feature csv
	prediction_df=make_test_df(bin_models, unknown_df, second_layer_size,max_angstrom,first_layer_size,staggering_value, features_to_include)
	if test==True:
		df_nom=RMSD_nom(unknown_df,second_layer_list_up,second_layer_list_down)
		df_nom=df_nom[["local_CA_nom","ID"]]
		#using third layer model we make final prediction
		unknown_df=unknown_df.merge(df_nom,on="ID")
	a=prediction_df.join(unknown_df,on=None)
	#print(a.columns)
	second_layer_list_up=np.arange(second_layer_size,max_angstrom,second_layer_size)
	second_layer_list_up=[str(x) for x in second_layer_list_up]
	second_layer_list_up.append(str(max_angstrom))
	#print(second_layer_list_up)
	if test==True:
		list_features=second_layer_list_up+features_to_include+["local_CA_nom"]
	else:
		list_features=second_layer_list_up+features_to_include

	a=a[list_features]
	a=a.drop(labels=["ID"], axis=1, errors='ignore')
	model= joblib.load(os.path.join(third_model_path))
	if test==True:
		nom=a[["local_CA_nom"]].to_numpy()
		a=a.drop(labels="local_CA_nom",axis=1).to_numpy()
	else:
		a.to_numpy()
	predictions=model.predict(a)

	if test==True:
		print("Predicted value: ",predictions,"Actual local C-alpha RMSD: ",nom)
		predictions=[predictions, nom]


	return (predictions)



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