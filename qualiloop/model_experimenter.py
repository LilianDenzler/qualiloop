#!/usr/bin/python3
"""/*************************************************************************
	Program:    Qualiloop Model Experimenter
	File:       modelmaker.py
	
	Version:    V1.1
	Date:       25.03.21
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
	The program contains all functions for an experiment run to build the best
	possible qualiloop model
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
import math 
import time
from datetime import datetime
#from pyinstrument import Profiler
from qualiloop import save_RMS_lib
from qualiloop import modelmaker
from qualiloop import visualizer
from qualiloop import preprocessor
from qualiloop import numerical_modelmaker
#import save_RMS_lib
#import modelmaker
#import visualizer

from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold
import xgboost as xgb


from sklearn.metrics import matthews_corrcoef
#import doctest
#import preprocessor
from sklearn.base import BaseEstimator
from sklearn.metrics import matthews_corrcoef

from math import sqrt
import argparse
import configparser
from datetime import datetime
from pathlib import Path
import warnings


def prepare_full_dataset_df(model_dir, log_dir, actual_dir, model_type, mode, params=None,full_dataset_path=None,save_prefix=None,run_dir=None):

	equal_params=None
	balanced_params=None
	staggered_params=None
	if model_type=="staggered":
		staggered_params=params
	if model_type=="equal":
		print("here")
		equal_params=params
	if model_type=="balanced":
		balanced_params=params

	#for convenience, open file, if full dataset has already been calculated
	if full_dataset_path:
		full_dataset_df=pd.read_csv(os.path.join(full_dataset_path),header=0, index_col=False)
		RMSD_num_df=full_dataset_df[mode]
		val_df=modelmaker.make_val_df(balanced_params, staggered_params,equal_params,RMSD_num_df,mode)
	#make a new full dataset (takes some time)
	else:
		full_dataset_df, val_df=modelmaker.get_full_dataset_df(model_dir, log_dir,actual_dir, mode, equal_params=equal_params, balanced_params=balanced_params,
										staggered_params=staggered_params,save_to_path=os.path.join(run_dir,"full_dataset.csv"))

	#####################
	#clean and prepare full dataset
	#clean dataset, i.e. delete rows with nans, delete columns with only one value
	#encode amino acids
	print(full_dataset_df)
	full_dataset_df=modelmaker.clean_dataset(full_dataset_df)
	#no_outlier_df=preprocessor.remove_outliers(full_dataset_df, run_dir)
	#scaling, etc
	full_dataset_df.to_csv(os.path.join(run_dir,"complete_df_not_encoded.csv"), index=None)
	print(full_dataset_df.columns)
	visualizer.feature_box_plot(full_dataset_df, "features_box_plot.png", run_dir)
	print("BEFORE SCALING",full_dataset_df)
	full_dataset_df=preprocessor.scaling(full_dataset_df,"robust", run_dir)
	print("AFTER SCALING",full_dataset_df)
	visualizer.feature_box_plot(full_dataset_df, "features_box_plot_after_standard_scaling.png", run_dir)
	full_dataset_df=clean_dataset(full_dataset_df, run_dir)
	
	full_dataset_df=encode_full_df(full_dataset_df, run_dir)
	full_dataset_df=clean_dataset(full_dataset_df, run_dir)

	for column in full_dataset_df.columns:
			try:
				full_dataset_df[column]=full_dataset_df[column].to_numpy().astype(np.float)
			except:
				del full_dataset_df[column]
				print("deleted the following column:")
	#set aside validation set
	complete_df=full_dataset_df
	#train=full_dataset_df.sample(frac=0.8,random_state=42)
	#validation_dataset_df=full_dataset_df.drop(train.index)
	try:
		del full_dataset_df["level_0"]
	except:
		pass
	return complete_df, val_df #validation_dataset_df, full_dataset_df, 


def clean_dataset(full_dataset_df, run_dir):
	full_dataset_df.replace("nan",None)
	full_dataset_df.dropna(inplace=True)

	mask=~full_dataset_df.isin([np.nan, np.inf, -np.inf, "nan", "None"])
	deleted_cols=[]
	for column in mask:
		mask_nans=mask.loc[mask[column] == False]
		if mask_nans.shape[0]>=60:                               #####################CHANGE TO SOMETHING HIGHER IF NOT DUMMY MODE
			del full_dataset_df[column]
			deleted_cols.append(column)
	with open(os.path.join(run_dir,"info.txt"), 'a') as file:
		file.write("Deleted the following columns as tey contained more than 60 Nan values: {}\n".format(deleted_cols))
	visualizer.data_preanalysis_nans(full_dataset_df, "dropped_cols_"+ str(deleted_cols), run_dir=run_dir)

	before_rows=full_dataset_df.shape[0]
	full_dataset_df=modelmaker.clean_dataset(full_dataset_df)
	after_rows=full_dataset_df.shape[0]
	print("Number of rows deleted in cleaning process: ", before_rows-after_rows)
	full_dataset_df = full_dataset_df[full_dataset_df.columns.drop(list(full_dataset_df.filter(regex='index')))]
	return full_dataset_df

def encode_full_df(full_dataset_df, run_dir):
	#encodings
	non_num_cols=[]
	for column in full_dataset_df.columns:
		if column=="ID":
			continue
		try:
			full_dataset_df[column]=full_dataset_df[column].astype(np.float)
		except:
			all_number_list = list(map(lambda x: str(x).isdigit(), full_dataset_df[column]))

			start = time.time()
			mixed=False
			while (time.time() - start < 3):
				for val in full_dataset_df[column]:
					if any(char.isdigit() for char in str(val))==True:
						mixed=True
						break

			if len(set(all_number_list)) != 1 or mixed==True:
				full_dataset_df=preprocessor.encode_tip(full_dataset_df,column)
				print(column, "is mixed")
				continue
			else:
				non_num_cols+=[column]
	non_num_cols=["tip_res"]                                                  #############in dummy mode!!!!!!!!!!!!!!!
	full_dataset_df=preprocessor.aa_encoding_all(full_dataset_df, non_num_cols, run_dir) #doesnt work
	
	return full_dataset_df



def baseline_models(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode,run_name=None,run_dir=None):
	#1.a no-opti models, rudimentary features, predict nom directly
	#delete all non-numerical columns to ensure nothing missed during encoding stage
	MCC_dic=modelmaker.make_models(pd.DataFrame(X_train.to_numpy()), pd.DataFrame(y_train.to_numpy()), pd.DataFrame(X_test.to_numpy()), pd.DataFrame(y_test.to_numpy()),models)
	models_torun=pd.DataFrame()
	if len(MCC_dic.values()) >=1:
		MCC_df=pd.DataFrame.from_dict(MCC_dic, orient = 'columns')
		max_value = max(MCC_df["MCC_mean"].to_list())  # maximum value
		max_keys=list(MCC_df.loc[(MCC_df['MCC_mean'] == max_value, "Model_name")])
		for model_name, model in models:
				if [model_name]==max_keys:
					models_torun["RMSD_nom"]=[nom_col]
					models_torun["model_name"]= max_keys
					models_torun["MCC"]= [max_value]
					out_model=model
	if run_name==None:
		run_name="baseline_nom"
	models_torun["descriptor"]=[run_name]
	with open(os.path.join(run_dir,"info.txt"), 'a') as file:
		file.write("Baseline model results: {}\n".format(MCC_df))
		file.write("Best of the baseline models : {}\n".format(models_torun))
	#models_viz(X_train, y_train.to_numpy(), X_test, y_test.to_numpy(), RMSD_mode, model, run_name,run_dir=run_dir)
	
	return models_torun, out_model

def baseline_models_opti(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode,run_name=None,run_dir=None):
	#1.b opti models, rudimentary features, predict nom directly
	#delete all non-numerical columns, as no encoding is done at this stage

	models_torun, MCC_df, MCC_df_opti, model=modelmaker.get_best_submodel_opti(X_train, X_test, y_train,y_test,models)
	models_torun["RMSD_nom"]=[nom_col]
	if run_name==None:
		run_name="baseline_opti_nom"
	models_torun["descriptor"]=[run_name]
	with open(os.path.join(run_dir,"info.txt"), 'a') as file:
		file.write("Baseline model optimized results: {}\n".format(MCC_df))
		file.write("Best of the baseline optimized models : {}\n".format(models_torun))
	#models_viz(X_train, y_train.to_numpy(), X_test, y_test.to_numpy(), RMSD_mode, model, run_name,run_dir=run_dir)

	return models_torun, model


def models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name, run_dir=None):
	
	visualizer.ROCAUC_curve(X_train, y_train, X_test, y_test, RMSD_mode, model,"ROCAUC.png", run_name=run_name,run_dir=run_dir)
	
	visualizer.precision_recall_curve(X_train, y_train, X_test, y_test, RMSD_mode, model,"precision_recall.png", run_name=run_name,run_dir=run_dir)
	
	visualizer.random_forest_features(X_train, y_train, RMSD_mode, "random_forest_features",run_name=run_name, run_dir=run_dir)
	
	visualizer.learning_curve(X_train, y_train, model,"learning_curve.png", run_name=run_name,run_dir=run_dir)
	visualizer.feature_importance(X_train, y_train,model,"feature_importance.png", run_name=run_name,run_dir=run_dir)
	visualizer.recursive_feature_elimination(X_train, y_train,model,"recursive_feature_elimination.png", run_name=run_name,run_dir=run_dir)


def dataset_preprocess_viz(full_dataset_df, RMSD_mode, bin_thresholds, val_df,RMSD_name,run_dir):
	#Visualize Class Distributions
	visualizer.class_distribution(full_dataset_df, RMSD_mode, bin_thresholds, RMSD_name+"_class_distributions.png",  run_dir=run_dir)
	visualizer.num_distribution(full_dataset_df, RMSD_mode,bin_thresholds,val_df, RMSD_name,RMSD_name+"_distribution.png", run_dir=run_dir)

	#Various visualizations for feature analysis and correlation analysis
	visualizer.single_feature_corr(full_dataset_df, RMSD_mode, RMSD_name+"_features_correlation.png",  run_dir=run_dir)
	try:
		visualizer.random_forest_features(full_dataset_df, RMSD_mode, RMSD_name+"_random_forest_rank.png",  run_dir=run_dir)
	except:
		pass
	try:
		visualizer.pearson_ranker(full_dataset_df, RMSD_mode, RMSD_name+"_pearson_ranker.png", run_dir=run_dir)
	except:
		pass
	try:
		visualizer.covariance_ranker(full_dataset_df, RMSD_mode, RMSD_name+"_covaraince_ranker.png",  run_dir=run_dir)
	except:
		pass
	try:
		visualizer.correlation_plot_all(full_dataset_df, RMSD_name+"_corr_plot.png",  run_dir=run_dir)
	except:
		pass
	try:
		visualizer.PCA_decomposition(full_dataset_df, RMSD_mode, RMSD_name+"_PCA_decomposition.png",  run_dir=run_dir)
	except:
		pass


def baseline_model_experimenter(X_train,y_train, X_test,y_test,RMSD_mode, del_y, models, models_torun_results, save_prefix,run_dir=None):
	#TEST PERFORMANCE OF BASELINE MODELS, RUDIMENTARY FEATURES
	X_train_rud=X_train[["length", "identity", "similarity"]]
	X_test_rud=X_test[["length", "identity", "similarity"]]
	nom_col= [col for col in RMSD_mode if "nom" in col][0]
	y_train=y_train[[nom_col]]
	y_test=y_test[[nom_col]]

	#no hyperparameter optimization
	models_torun_baseline_rud_row, model=baseline_models(X_train_rud, X_test_rud, y_train, y_test, nom_col, models, RMSD_mode,run_name="baseline_rud",  run_dir=run_dir)
	models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
	print(models_torun_results)
	models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_rud_non_opti"+".csv"))
	#models_viz(X_train_rud, y_train, X_test_rud, y_test, RMSD_mode, model, run_name="baseline_rud" ,run_dir=run_dir)

	rf_prediction=model.predict(X_test_rud.to_numpy())
	print("PREDICTION",rf_prediction)
	rf_prediction_str=pd.DataFrame(rf_prediction)
	rf_prediction_str=pd.concat([rf_prediction_str, y_test], axis=1)
	rf_prediction_str.to_csv(os.path.join(run_dir,"Predictions_rud_non_opti"+".csv"))
	MCC_dic=pd.DataFrame([matthews_corrcoef(y_test,rf_prediction, sample_weight=None )])
	#MCC_dic=modelmaker.make_models(X_train_rud, y_train, X_test_rud, y_test,[(type(model).__name__, model)],only_score=True)
	#MCC_dic=pd.DataFrame.from_dict(MCC_dic)
	MCC_dic.to_csv(os.path.join(run_dir,"Validation_rud_non_opti"+".csv"))

	#with hyperparameter optimization
	models_torun_baseline_rud_row, model=baseline_models_opti(X_train_rud, X_test_rud, y_train, y_test, nom_col, models, RMSD_mode,run_name="baseline_rud_opti",run_dir=run_dir)
	models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
	print(models_torun_results)
	models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_rud_opti"+".csv"))
	#models_viz(X_train_rud, y_train, X_test_rud, y_test, RMSD_mode, model, run_name="baseline_rud_opti",run_dir=run_dir)
	full_y=np.concatenate((np.unique(y_train), np.unique(y_test)))
	num_class=len(np.unique(full_y))
	"""if type(model).__name__=="XGBClassifier":
		print("XGB!!!!!!!")
		param={}
		param['objective']='multi:softmax'
		param['num_class']=num_class+1
		model=model.set_params(**param)
		dMatrixTest = xgb.DMatrix(X_test_rud, label=y_test)
		rf_prediction=model.predict(X_test_rud)
	else:
		rf_prediction=model.predict(X_test_rud.to_numpy())"""
	rf_prediction=model.predict(X_test_rud)
	print("PREDICTION",rf_prediction)
	rf_prediction_str=pd.DataFrame(rf_prediction)
	rf_prediction_str=pd.concat([rf_prediction_str, y_test], axis=1)
	rf_prediction_str.to_csv(os.path.join(run_dir,"Predictions_rud_opti"+".csv"))
	MCC_dic=pd.DataFrame([matthews_corrcoef(y_test,rf_prediction, sample_weight=None )])
	#MCC_dic=modelmaker.make_models(X_train_rud, y_train, X_test_rud, y_test,[(type(model).__name__,model)],only_score=True)
	#MCC_dic=pd.DataFrame.from_dict(MCC_dic)
	MCC_dic.to_csv(os.path.join(run_dir,"Validation_rud_opti"+".csv"))

	#TEST PERFORMANCE OF BASELINE MODELS, ALL FEATURES
	#no hyperparameter optimization
	models_torun_baseline_rud_row, model=baseline_models(X_train, X_test, y_train, y_test, nom_col, models,RMSD_mode,run_name="baseline_all",run_dir=run_dir)
	models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
	print(models_torun_results)
	models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_all_non_opti"+".csv"))
	#models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name="baseline_all",run_dir=run_dir)

	rf_prediction=model.predict(X_test.to_numpy())
	print("PREDICTION",rf_prediction)
	rf_prediction_str=pd.DataFrame(rf_prediction)
	rf_prediction_str=pd.concat([rf_prediction_str, y_test], axis=1)
	rf_prediction_str.to_csv(os.path.join(run_dir,"Predictions_all_non_opti"+".csv"))
	MCC_dic=pd.DataFrame([matthews_corrcoef(y_test,rf_prediction, sample_weight=None )])
	#MCC_dic=modelmaker.make_models(X_train_rud, y_train, X_test_rud, y_test,[(type(model).__name__,model)],only_score=True)
	#MCC_dic=pd.DataFrame.from_dict(MCC_dic)
	MCC_dic.to_csv(os.path.join(run_dir,"Validation_all_non_opti"+".csv"))


	#with hyperparameter optimization
	models_torun_baseline_rud_row, model=baseline_models_opti(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode, run_name="baseline_all_opti",run_dir=run_dir)
	models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
	print(models_torun_results)
	models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_all_opti"+".csv"))
	#models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name="baseline_all_opti",run_dir=run_dir)

	rf_prediction=model.predict(X_test)
	print("PREDICTION",rf_prediction)
	rf_prediction_str=pd.DataFrame(rf_prediction)
	rf_prediction_str=pd.concat([rf_prediction_str, y_test], axis=1)
	rf_prediction_str.to_csv(os.path.join(run_dir,"Predictions_all_opti"+".csv"))
	MCC_dic=pd.DataFrame([matthews_corrcoef(y_test,rf_prediction, sample_weight=None )])
	#MCC_dic=modelmaker.make_models(X_train_rud, y_train, X_test_rud, y_test,[(type(model).__name__,model)],only_score=True)
	#MCC_dic=pd.DataFrame.from_dict(MCC_dic)
	MCC_dic.to_csv(os.path.join(run_dir,"Validation_all_opti"+".csv"))


	return models_torun_results

def num_predict(X_train,y_train, X_test,y_test,RMSD_name, val_df_mode, models_torun_results, save_prefix,run_dir=None, already_trained=False):
	#https://github.com/WillKoehrsen/Data-Analysis/blob/master/prediction-intervals/prediction_intervals.ipynb
	models_torun_results=pd.DataFrame()
	nom_col=[col for col in y_train.columns if "nom" in col and RMSD_name in col]
	y_train_nom=y_train[nom_col]
	y_test_nom=y_test[nom_col]
	y_train=y_train[[RMSD_name]]
	y_test=y_test[RMSD_name]
	print("gradientBoost_quantile")
	predictions, within_df, model_MCCs=numerical_modelmaker.gradientBoost_quantile(X_train, y_train, X_test,y_test, val_df_mode,RMSD_name,run_dir=run_dir, already_trained=already_trained)
	lower_MCC=model_MCCs.loc["lower_model"]["MCC"]
	mid_MCC=model_MCCs.loc["mid_model"]["MCC"]
	upper_MCC=model_MCCs.loc["upper_model"]["MCC"]
	models_torun_results["lower_model"]=[lower_MCC]
	models_torun_results["mid_model"]=[mid_MCC]
	models_torun_results["upper_model"]=[upper_MCC]

	print("linear_regression_model")
	#,MCC, MCC_train
	linear_regression_preds,linear_regression_description, test_mse, test_mae, test_rmse=numerical_modelmaker.linear_regression_model(X_train, y_train,
																																		X_test,y_test, val_df_mode, run_dir=run_dir, already_trained=already_trained)
	print(linear_regression_description)
	print(test_mse)
	print(test_mae)
	print(test_rmse)
	#models_torun_results["linear_regression_model"]=[MCC]

	MCC, MCC_train=numerical_modelmaker.keras_model(X_train, y_train,X_test,y_test,"keras_model.png",val_df_mode,run_dir=run_dir, already_trained=already_trained)
	models_torun_results["keras_model"]=MCC
	MCC, MCC_train=numerical_modelmaker.quantile_loss_keras(X_train, y_train,X_test,y_test, "quantile_loss_NN",val_df_mode, run_dir=run_dir, already_trained=already_trained)
	models_torun_results["quantile_loss_keras_model"]=[MCC]
	print("XGBOOSTREGRESSION")
	wrong_count_test, wrong_count_train, MCC, MCC_train=numerical_modelmaker.xgboost_model(X_train, y_train,X_test,y_test, val_df_mode, 
																			"xgboost",run_dir=run_dir, already_trained=already_trained)
	
	models_torun_results["XGBRegressor_model"]=[MCC]

	print(models_torun_results)

	return models_torun_results


def full_experimentation(model_dir, log_dir, actual_dir, model_type, mode, params=None, full_dataset_path=None, manual_sele=None, save_prefix=None,run_dir=None):
	
	#prepare full dataset with chosen target-y and all features. Also extract validation set, to be used for final testing. 
	if run_dir==None:
		run_dir=datetime.today().strftime('%Y_%m_%d_%H-%M-%S_qualiloop_run')
	if save_prefix==None:
		save_prefix=""
	if not os.path.exists(os.path.join(run_dir)):
			os.makedirs(os.path.join(run_dir))

	with open(os.path.join(run_dir,"info.txt"), 'a') as file:
		file.write("""START OF RUN\nmodel_type: {},\nmode: {},\nparams: {},\nfull_dataset_path: {},\nsave_prefix: {},\nrun_dir: {}\n""".format(model_type, mode, params, full_dataset_path,save_prefix,run_dir))

	complete_df, val_df=prepare_full_dataset_df(model_dir, log_dir, actual_dir,model_type, mode, params=params, full_dataset_path=full_dataset_path, 
																											save_prefix=save_prefix, run_dir=run_dir)
	
	try:
		complete_df.rename(columns = {'local_CA_bin5.199999999999999':'local_CA_bin5.2'}, inplace = True)
	except:
		pass
	#define the model types used 
	models=modelmaker.all_classifier_types(model_list=None)
	run_dir_all=run_dir
	#start of model-making process
	for RMSD_name in mode:
		print(RMSD_name)
		print(val_df["bin_thresholds"].to_list())
		feature_cols=[col for col in complete_df.columns if "nom" not in col and "bin" not in col and "local" not in col and "global" not in col]
		print("feature_cols",feature_cols)
		RMSD_y_cols_nombin=val_df["bin_thresholds"].to_list()[0]
		list_append=[]
		for x in RMSD_y_cols_nombin:
			list_append.append(RMSD_name+"_bin"+str(x))
		RMSD_y_cols_nombin=list_append
		RMSD_y_cols_nombin.append(RMSD_name+"_nom")
		RMSD_y_cols=RMSD_y_cols_nombin.copy()
		RMSD_y_cols.append(RMSD_name)
		print(RMSD_y_cols)

		#TEST function to ensure rounding was done on full_dataset
		if RMSD_y_cols==[col for col in complete_df.columns if col not in feature_cols] and RMSD_y_cols_nombin==[ele for ele in RMSD_y_cols if ele not in mode]:
			pass
		else:
			pass
			#raise Warning("The dataset contains unrounded or differing values for the binary cut-off values.")

		print("zguhkJ",RMSD_y_cols_nombin)
		run_dir=os.path.join(run_dir_all,RMSD_name)
		if not os.path.exists(os.path.join(run_dir,"graphs")):
			os.makedirs(os.path.join(run_dir,"graphs"))
		#set aside validation set
		print(complete_df.columns)
		train,vali = train_test_split(complete_df, stratify=complete_df[RMSD_name+"_nom"], test_size=0.2)
		print(train.columns)
		X_train=train[feature_cols]
		y_train=train[RMSD_y_cols]
		X_test=vali[feature_cols]
		y_test=vali[RMSD_y_cols]
		print("RMSD_y_cols", RMSD_y_cols)
		print(X_train.shape, y_train.shape,"correct?")
		val_df_mode=val_df.loc[RMSD_name , : ]
		#dataset_preprocess_viz(complete_df, RMSD_y_cols_nombin, val_df_mode["bin_thresholds"], val_df_mode,RMSD_name,run_dir=run_dir)

		#feature selection
		with open(os.path.join(run_dir,"info.txt"), 'a') as file:
			file.write("FEATURE SELECTION\n")
			file.write("Feature selection with max number of features : {}\n".format(30))
		print(X_train,y_train, RMSD_y_cols_nombin)
		selection_dic=preprocessor.feature_selection(X_train,y_train, RMSD_y_cols_nombin, 30, run_dir, manual_sele=manual_sele)
		run_dir_og=run_dir
		X_train_og=X_train
		X_test_og=X_test
		for selection in selection_dic.keys():
			run_dir=run_dir_og
			X_train=X_train_og
			X_test=X_test_og
			run_dir=os.path.join(run_dir,str(selection))
			if not os.path.exists(os.path.join(run_dir)):
				os.makedirs(os.path.join(run_dir))
			all_cols=np.concatenate([selection_dic[selection]],axis=0)
			sele=pd.DataFrame()
			sele["sele"]=selection_dic[selection]
			sele.to_csv(os.path.join(run_dir,"selected_features.csv"))
			X_test=X_test[selection_dic[selection]]
			y_test=y_test[RMSD_y_cols]
			X_train=X_train[selection_dic[selection]]
			y_train=y_train[RMSD_y_cols]
			if not os.path.isdir(run_dir):
				os.makedirs(run_dir)
			with open(os.path.join(run_dir,"info.txt"), 'a') as file:
				file.write("STARTING RUN WITH SELECTION\n".format(selection))
			with open(os.path.join(run_dir,"selected_features.txt"), 'w') as file:
				file.writelines(selection_dic[selection])
			vif=preprocessor.calc_vif_score(X_train[selection_dic[selection]], RMSD_y_cols)
			vif.to_csv(os.path.join(run_dir,"vif_features_common_sele.csv"))
			with open(os.path.join(run_dir,"info.txt"), 'a') as file:
				file.write("STARTING MODEL GENERATION FOR RMSD MODE: {}\n".format(RMSD_y_cols_nombin))

			#save to this subfolder to keep information together
			val_df_mode.to_csv(os.path.join(run_dir,"val_df.csv"))
			vali.to_csv(os.path.join(run_dir,"validation_dataset.csv"), index=None)
			train.to_csv(os.path.join(run_dir,"train_dataset.csv"), index=None)
			vali.to_csv(os.path.join(run_dir,"test_dataset.csv"), index=None)
			complete_df.to_csv(os.path.join(run_dir,"complete_dataset.csv"), index=None)
			#make dataset with only the current RMSD_y_cols_nombin's target y columns
			X_train=modelmaker.clean_dataset(X_train)
			X_test=modelmaker.clean_dataset(X_test)

		  
			#MODEL MAKER FULL RUN#######################################################################################################
			models_torun_results=pd.DataFrame() #dataframe containing final MCC of the best model for each specified run

			####################
			#NUMERICAL PREDICTOR
			####################
			#models_torun_results_num=num_predict(X_train,y_train, X_test,y_test,RMSD_name, val_df_mode,models_torun_results, save_prefix,run_dir=run_dir, already_trained=False)
			#print(models_torun_results_num)
			#quit()
			################
			#BASELINE MODELS
			################ 
			#models_torun_results=baseline_model_experimenter(X_train,y_train, X_test,y_test,RMSD_y_cols_nombin, RMSD_y_cols, models, models_torun_results,save_prefix,  run_dir=run_dir)
				
			################
			#STACKED MODELS
			################
			with open(os.path.join(run_dir,"info.txt"), 'a') as file:
				file.write("STARTING STACKED MODEL\n")
			print(run_dir)
			(y_train_staggered,y_test_staggered,models_torun_double_staggered_all, model_all, 
				models_torun_only_sums_double_staggered, model_only_sums,MCC_of_bins,X_train_with_preds,
				X_test_with_preds,second_layer_df_train,second_layer_df_test)= modelmaker.double_staggered_run(X_train,y_train, RMSD_y_cols,models, 
																																			RMSD_y_cols_nombin, val_df_mode, [model_type], run_dir=run_dir)
			models_torun_results=pd.concat([models_torun_results, models_torun_double_staggered_all, MCC_of_bins])
			models_torun_results=pd.concat([models_torun_results, models_torun_only_sums_double_staggered, MCC_of_bins])
			models_torun_results.to_csv(os.path.join(run_dir,"OVERNIGHT_RESULTS_stacked"+".csv"))
			print(X_train_with_preds,y_train_staggered)
			print(X_test_with_preds,y_test_staggered)
			#models_viz(X_train_with_preds, y_train_staggered, X_test_with_preds, y_test_staggered, RMSD_y_cols_nombin, model_all, "stacked_models_all", run_dir=run_dir)
			models_viz(second_layer_df_train, y_train_staggered,second_layer_df_test, y_test_staggered, RMSD_y_cols_nombin, model_only_sums, "stacked_models_only_sums", run_dir=run_dir)



			try:
				scores = model.score(X_test.to_numpy(), y_test.to_numpy)
				scores_no_sums = model_only_sums.score(X_test.to_numpy(), y_test.to_numpy())
				print(scores)
				print(scores_only_sums)
				with open(os.path.join(os.path.join(run_dir,"SCORES_stacked"+".csv")), 'a') as file:
					file.write(str(scores))
					file.write(str(scores_only_sums))
			except:
				pass

			run_dir=os.path.join(run_dir+"_no_opti")
			models_torun_results=pd.DataFrame()
			(y_train_staggered,y_test_staggered,models_torun_double_staggered_all_no_opti, 
				model_all_no_opti, models_torun_double_staggered_only_sums_no_opti, model_only_sums_no_opti,
				MCC_of_bins,X_train_with_preds,X_test_with_preds,second_layer_df_train,second_layer_df_test)=modelmaker.double_staggered_run(X_train,y_train, 
																																			RMSD_y_cols,models, RMSD_y_cols_nombin, val_df_mode, 
																																			model_type, run_dir=run_dir, no_opti=True)
			models_torun_results=pd.concat([models_torun_results, models_torun_double_staggered_no_opti, MCC_of_bins])
			models_torun_results=pd.concat([models_torun_results, models_torun_double_staggered_only_sums_no_opti, MCC_of_bins])
			models_torun_results.to_csv(os.path.join(run_dir,"OVERNIGHT_RESULTS_stacked"+".csv"))
			#models_viz(X_train_with_preds, y_train_staggered, X_test_with_preds, y_test_staggered, RMSD_y_cols_nombin, model_all, "stacked_models_all", run_dir=run_dir)
			models_viz(second_layer_df_train, y_train_staggered,second_layer_df_test, y_test_staggered, RMSD_y_cols_nombin, model_only_sums, "stacked_models_only_sums", run_dir=run_dir)

				
			#MCC_dic_only_sums=modelmaker.make_models(X_train, y_train, X_test, y_test,[(type(model_only_sums).__name__,model_only_sums)],only_score=True)
			#MCC_dic_only_sums=pd.DataFrame.from_dict(MCC_dic_only_sums)
			#MCC_dic_only_sums.to_csv(os.path.join(run_dir,"Validation_stacked_only_sums"+".csv"))
			#MCC_dic=modelmaker.make_models(X_train, y_train, X_test, y_test,[(type(model).__name__,model)],only_score=True)
			#MCC_dic=pd.DataFrame.from_dict(MCC_dic)
			#MCC_dic.to_csv(os.path.join(run_dir,"Validation_stacked"+".csv"))

			#test on validation set
def check_args(parser,settings):
	for key,val in settings.items():
		try:
			val=eval(val)
			settings.update({key:val})
		except:
			pass
	if settings["type"]=='equal':
		if settings["ep"]==None:
			parser.error("-equal requires -equal_params (-ep)")
		else:
			params=settings["ep"]
			try:
				params=np.array(params,float)
				if type(params[0])!=np.array:
					params=[params]
			except:
				parser.error("-ep input must be numerical")

	elif settings["type"]=='staggered':
		if settings["sp"]==None:
			parser.error("-staggered requires -staggered_params (-sp)")
		else:
			params=settings["sp"]
			try:
				params=np.array(params,float)
				if type(params[0])!=np.array:
					params=[params]
			except:
				parser.error("-sp input must be numerical")
	elif settings["type"]=='balanced':
		if settings["bp"]==None:
			parser.error("-balanced requires -balanced_params (-bp)")
		else:
			if type(params[0])!=np.array and type(params[0])!=list:
				params=[params]
			params=settings["bp"]
	if settings["type"]==None:
		parser.error("the following arguments is required: type")
	if settings["mode"]==None:
		parser.error("the following arguments is required: mode")
	if type(settings["mode"])==str:
		settings["mode"]=settings["mode"].split(",")
	elif type(settings["mode"])==list:
		if len(settings["mode"])>4:
			parser.error("Too many modes. Maximum of four: local_CA, local_AA, global_CA, global_AA.")
		settings["mode"] = [item for sublist in settings["mode"] for item in sublist]
	print("HEIHE", params)
	print("nums",len(params),len(settings["mode"]))

	print("params",params)
	if len(params)==1:
		if len(params)==len(settings["mode"]):
			pass
		else:
			params=params*len(settings["mode"])
			warnings.warn("Only one set of parameters given for multiple modes. The same parameters are used for each mode.")
	else:
		try:
			params=[params]
			if len(params)==1:
				params=params*len(settings["mode"])
				warnings.warn("Only one set of parameters given for multiple modes. The same parameters are used for each mode.")

		except:
			parser.error("Unequal number of parameter sets and modes. Please give parameters for each mode or one set of parameters for all modes.")
	print(settings["mode"])

	try:
		if settings["full"]=="":
			settings["full"]=None
	except:
		settings["full"]=None
		
	return settings,params

def main():

	conf_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,add_help=False)
	conf_parser.add_argument("-config", "--c",dest="config", help="Specify config file", metavar="FILE")
	args, remaining_argv = conf_parser.parse_known_args()

	settings = {}

	if args.config:
		config = configparser.SafeConfigParser()
		config.read([args.config])
		print(config.items("model_1"))
		settings.update(dict(config.items("model_1")))

	parser = argparse.ArgumentParser(epilog="For any comments or inquiries please contact zcbtlm0@ucl.ac.uk",#fromfile_prefix_chars='@', 
		formatter_class=argparse.RawDescriptionHelpFormatter,parents=[conf_parser],#, argparse.ArgumentDefaultsHelpFormatter
		description=('''\
			Please do not mess up this text!
			--------------------------------
					I have indented it
					exactly the way
					I want it
			''')) #textwrap.dedent
	parser.set_defaults(**settings)
	verbosity = parser.add_mutually_exclusive_group()
	verbosity.add_argument("-v", "--verbose", action="store_true")
	verbosity.add_argument("-q", "--quiet", action="store_true")

	parser.add_argument('m', nargs='?',type=Path, metavar="model_dir_path", help='path of directory where all abYmod-generated files are stored')
	parser.add_argument('l', nargs='?', type=Path,metavar="log_dir_path", help='path of directory where all abYmod-generated log-files are stored')
	parser.add_argument('a', nargs='?', type=Path,metavar="actual_dir_path", help='path of directory where all actual PDB files are stored')

	parser.add_argument('-type' , metavar='type',choices=['staggered', 'equal', 'balanced'], help='',default=None)
	parser.add_argument('-mode',  nargs='+', action='append',type=str,metavar="mode", help='', default=None)
	
	parser.add_argument('-sp', nargs='+', action='append',metavar="vals", help='upper threshold, lower threshold and number of binary classes for global AA ')
	parser.add_argument('-ep',nargs='+', action='append',metavar="thresholds", help='upper threshold, lower threshold and number of binary classes for local AA ')
	parser.add_argument('-bp', nargs='+',action='append',type=int, metavar="vals", help='upper threshold, lower threshold and number of binary classes for global AA ')
	
	parser.add_argument('-full', default=None, type=Path,nargs='+',  help='finished dataset path')
	parser.add_argument('-name',default=None, type=Path,help='save the project run under this name.')
	parser.add_argument('-r',  default=None,type=Path,help='save your project in this directory.')
	parser.add_argument('-sele', default=None,nargs='+',help='names of manually selected features')

	#parser.add_argument('-config',metavar="FILE", dest="config_path", help='path to config file.')
	args2 = vars(parser.parse_args())
	print(args)
	
	settings.update({k: v for k, v in args2.items() if v is not None})
	settings.update({k: v for k, v in args2.items() if v is None and k not in settings.keys()}) 
	
	settings,params=check_args(parser,settings)
	full_experimentation(settings["m"],settings["l"], settings["a"],settings["type"], settings["mode"], params=params,full_dataset_path=settings["full"],
		manual_sele=settings["sele"], save_prefix=settings["name"],run_dir=settings["r"])


def full_setup():
	conf_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,add_help=False)
	conf_parser.add_argument("-config", "--c",dest="config", help="Specify config file", metavar="FILE")
	args, remaining_argv = conf_parser.parse_known_args()

	settings = {}
	config = configparser.SafeConfigParser()
	config.read([args.config])
	print(args.config)
	for section in config.sections():
		config_dic=dict(config.items(section))
		settings,params=check_args(conf_parser,config_dic)
		print(settings)
		print(params)
		full_experimentation(settings["m"],settings["l"], settings["a"],settings["type"], settings["mode"], params=params,full_dataset_path=settings["full"],
		manual_sele=settings["sele"], save_prefix=settings["name"],run_dir=settings["r"])


if __name__ == '__main__':
	#main()
	full_setup()
	"""model_dir=sys.argv[1]
	log_dir=sys.argv[2]
	actual_dir=sys.argv[3]
	equal_params=[[1,7,8], [1,8,6],[2,10,8],[1,10,8]] #"local_AA", "local_CA", "global_AA", "global_CA"
	bins_in_class=[2,2,2,2]
	nr_of_balanced_bins=10

	manual_sele=[['tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff'], 
	["similarity","length","identity"],
	['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', "contacts_all"], 
	['blosum_dist', 'contacts_all', 'Screlacc', 'Access', 'simlength', 'Scacc',  'identity', 'Relacc', 'tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff'], 
	['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', 'simlength', "PCAtip_res_blosum_621_blosum_62_","PCAtip_res_blosum_622_blosum_62_","PCAtip_res_blosum_623_blosum_62_"],
	['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', 'simlength', "PCAtip_res_NLF1_NLF_col","PCAtip_res_NLF2_NLF_col","PCAtip_res_NLF3_NLF_col"],
	['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', 'simlength', "res_charge_physiochem_4d_","res_sc_nr_physiochem_4d_","res_compactness_physiochem_4d_"]]


	#full_dataset_path="~/sync_project/WWW/CDRH3loop/CDRH3loop/other_stuff/staggered_try_out/staggered_try_out.csv"
	full_dataset_path="~/sync_project/WWW/CDRH3loop/staggered_0.2_5_1.5_2_new_copy/staggered_0.2_5_1.5_2_new.csv"
	bin_thresholds=[1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9]
	nom_values=[1,2,3]
	threshlow=[0,2,4]
	threshup=[2,4,100]
	values=[1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9]
	val_df=pd.DataFrame()
	for RMSD_mode in ["local_CA","local_AA","global_CA","global_AA"]:
		row = pd.Series({"nom_label":nom_values, "lower_threshold":threshlow, "upper_threshold":threshup, "bin_thresholds":values},name=RMSD_mode)
		val_df = val_df.append(row)

	full_dataset_pack=[full_dataset_path, bin_thresholds, val_df]
	save_prefix="staggered_0.2_5_1.5_2_new_copy"
	staggered_params=[[0.2,6,1.5,2],[0.2,5,1.5,2],[0.2,9,1.5,2],[0.2,9,1.5,2]] #[staggering_value, max_angstrom, first_layer_size, second_layer_size]
	run_dir=os.path.join('./', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, staggered=True, balanced=False,staggered_params=staggered_params,
								full_dataset_path=full_dataset_pack, manual_sele=manual_sele,save_prefix=save_prefix, run_dir=run_dir, fast_mode=None)
	"""
	"""
	full_dataset_path="/serv/www/html_lilian/CDRH3loop/CDRH3loop/staggered_0.1_5_1_2_global_CA/staggered_0.1_5_1_2_global_CA.csv"
	bin_thresholds=[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,
	3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,
	7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9]
	val_df={"nom_label":[1.0,2.0,3.0], "lower_threshold":[0,2,4,6,8], "upper_threshold":[2,4,6,8,100]}
	full_dataset_pack=[full_dataset_path, bin_thresholds, val_df]
	
	save_prefix="staggered_0.1_5_1_2_global_CA"
	staggered_params=[[0.1,6,1.5,2],[0.1,5,1,2],[0.1,9,1,2],[0.1,9,1,2]] #[staggering_value, max_angstrom, first_layer_size, second_layer_size]
	run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, staggered=True,staggered_params=staggered_params,
								full_dataset_path=full_dataset_pack, manual_sele=manual_sele,save_prefix=save_prefix, run_dir=run_dir)
								
	full_dataset_pack=None
	save_prefix="staggered_0.1_5_1.5_1"
	staggered_params=[[0.1,6,1.5,1],[0.1,5,1.5,1],[0.1,9,1.5,1],[0.1,9,1.5,1]] #[staggering_value, max_angstrom, first_layer_size, second_layer_size]
	run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, staggered=True,staggered_params=staggered_params,
								full_dataset_path=None, manual_sele=manual_sele,save_prefix=save_prefix, run_dir=run_dir)
	"""

	full_dataset_path=None
	equal_params=[[0,7,14], [0,8,16],[2,10,16],[2,10,16]] #"local_AA", "local_CA", "global_AA", "global_CA"
	bins_in_class=[2,2,2,2]
	nr_of_balanced_bins=[10,10,10,10]
	save_prefix="0816_2_uniform"
	run_dir=os.path.join('./', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, balanced=False,
								full_dataset_path=full_dataset_path, manual_sele=manual_sele, save_prefix=save_prefix, run_dir=run_dir)

	full_dataset_path=None
	equal_params=[[1,7,9], [1,8,12],[2,10,18],[1,10,9]] #"local_AA", "local_CA", "global_AA", "global_CA"
	bins_in_class=[3,3,3,3]
	nr_of_balanced_bins=[10,10,10,10]
	save_prefix="1812_3_uniform"
 
	run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, balanced=False,
								full_dataset_path=full_dataset_path, save_prefix=save_prefix, run_dir=run_dir)

	"""full_dataset_path=None
	equal_params=[[1,7,9], [1,8,12],[2,10,9],[1,10,9]] #"local_AA", "local_CA", "global_AA", "global_CA"
	bins_in_class=[4,4,4,4]
	nr_of_balanced_bins=[10,10,10,10]
	save_prefix="1812_4_uniform"
	run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, balanced=False,
								full_dataset_path=full_dataset_path, save_prefix=save_prefix, run_dir=run_dir)

	full_dataset_path=None
	equal_params=[[1,7,8], [1,8,6],[2,10,8],[1,10,8]] #"local_AA", "local_CA", "global_AA", "global_CA"
	bins_in_class=[2,2,2,2]
	nr_of_balanced_bins=[10,10,10,10]
	save_prefix="186_2_balanced"
	#now = datetime.now()
	#save_prefix=now.strftime("%d%b_%H%M%S")+save_prefix
	run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
	if not os.path.isdir(run_dir):
		os.makedirs(run_dir)
	full_experimentation(model_dir, log_dir, actual_dir, equal_params,bins_in_class,nr_of_balanced_bins, balanced=True,
								full_dataset_path=full_dataset_path, save_prefix=save_prefix, run_dir=run_dir)
	"""