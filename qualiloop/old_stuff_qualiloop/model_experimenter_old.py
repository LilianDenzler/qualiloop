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
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./')
sys.path.append('~/sync_project/WWW/CDRH3loop')
import os
os.environ['MODULE'] = '/serv/www/html_lilian/CDRH3loop'
module_path = os.getenv('MODULE')
import pandas as pd
import numpy as np
import joblib
import math 
import time
from datetime import datetime
from pyinstrument import Profiler

import save_RMS_lib
import modelmaker
import visualizer

from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold

import doctest
import preprocessor

def prepare_full_dataset_df(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins,full_dataset_path=None, fast_mode=None, save_prefix=None, balanced=False, staggered=False, stagger_params_all=None, run_dir=None):
   #for convenience, open file, if full dataset has already been calculated
   if full_dataset_path:
      full_dataset_df=pd.read_csv(os.path.join(full_dataset_path[0]),header=0, index_col=False)
      bin_thresholds= full_dataset_path[1]
      val_dic= full_dataset_path[2]
   #make a new full dataset (takes some time)
   else:
      full_dataset_df, bin_thresholds, val_dic=modelmaker.get_full_dataset_df(model_dir, log_dir,actual_dir,threshold_bins_staggered=threshold_bins_staggered, 
         step_size_nom_staggered=step_size_nom_staggered,save_to_path=os.path.join(run_dir,save_prefix+".csv"),
         nr_of_balanced_bins=nr_of_balanced_bins, balanced=balanced, fast_mode=fast_mode, staggered=staggered, stagger_params_all=stagger_params_all)

   #####################
   #clean and prepare full dataset
   #clean dataset, i.e. delete rows with nans, delete columns with only one value
   #encode amino acids
   full_dataset_df=modelmaker.clean_dataset(full_dataset_df)
   no_outlier_df=preprocessor.remove_outliers(full_dataset_df, run_dir)
   #scaling, etc
   visualizer.feature_box_plot(full_dataset_df, "features_box_plot.png", run_dir)
   print("BEFORE SCALING",full_dataset_df)
   full_dataset_df=preprocessor.scaling(full_dataset_df,"robust", run_dir)
   print("AFTER SCALING",full_dataset_df)
   visualizer.feature_box_plot(full_dataset_df, "features_box_plot_after_standard_scaling.png", run_dir)
   full_dataset_df=clean_dataset(full_dataset_df, run_dir)
   
   full_dataset_df=encode_full_df(full_dataset_df, run_dir)
   for column in full_dataset_df.columns:
         try:
            full_dataset_df[column]=full_dataset_df[column].to_numpy().astype(np.float)
         except:
            del full_dataset_df[column]
            print("deleted the following column:")
            #input(column)
   #set aside validation set
   full_dataset_df.to_csv(os.path.join(run_dir,save_prefix+"_encoded.csv"), index=None)
   complete_df=full_dataset_df
   train=full_dataset_df.sample(frac=0.8,random_state=42)
   validation_dataset_df=full_dataset_df.drop(train.index)
   try:
      del full_dataset_df["level_0"]
   except:
      pass
   return validation_dataset_df, full_dataset_df, complete_df, bin_thresholds, val_dic

def encode_provided():
   full_dataset_df=pd.read_csv("/serv/www/html_lilian/CDRH3loop/CDRH3loop/qualiloop/uniform_encoded.csv")

   #for col in full_dataset_df.columns:
   #   if "nom" in col:
   #      full_dataset_df[col] = full_dataset_df[col].subtract(1)

   bin_thresholds=[1.0,1.7142857142857144,2.428571428571429,3.142857142857143,3.857142857142857,4.571428571428571,5.285714285714286,6.0]
   val_dic={"nom_label":[1.0,2.0,3.0,4.0, 5.0], "lower_threshold":[0,1.7142857142857144,3.142857142857143,4.571428571428571,6], "upper_threshold":[1.7142857142857144,3.142857142857143,4.571428571428571,6,100]}
   
   complete_df=full_dataset_df
   train=full_dataset_df.sample(frac=0.8,random_state=42)
   validation_dataset_df=full_dataset_df.drop(train.index)
   try:
      del full_dataset_df["level_0"]
   except:
      pass
   return validation_dataset_df, full_dataset_df, complete_df, bin_thresholds, val_dic

def clean_dataset(full_dataset_df, run_dir):
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
   MCC_dic=modelmaker.make_models(X_train, y_train, X_test, y_test,models)
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
   models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name,run_dir=run_dir)
   
   return models_torun, out_model

def baseline_models_opti(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode,run_name=None,run_dir=None):
   #1.b opti models, rudimentary features, predict nom directly
   #delete all non-numerical columns, as no encoding is done at this stage

   models_torun, MCC_df, MCC_df_opti, model=modelmaker.get_best_submodel_opti(X_train, X_test, y_train, y_test,models)
   models_torun["RMSD_nom"]=[nom_col]
   if run_name==None:
      run_name="baseline_opti_nom"
   models_torun["descriptor"]=[run_name]
   with open(os.path.join(run_dir,"info.txt"), 'a') as file:
      file.write("Baseline model optimized results: {}\n".format(MCC_df))
      file.write("Best of the baseline optimized models : {}\n".format(models_torun))
   models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name,run_dir=run_dir)

   return models_torun, model

def prepare_for_train(full_dataset_df, RMSD_mode, del_y,included_features=None):
   for column in full_dataset_df.columns:
         try:
            full_dataset_df[column]=full_dataset_df[column].to_numpy().astype(np.float)
         except:
            del full_dataset_df[column]
            print("deleted the following column:")
            #input(column)

   nom_col= [col for col in RMSD_mode if "nom" in col][0]
   if included_features!=None:
      full_dataset_df=full_dataset_df[included_features+[nom_col]] 
   elif included_features==None:
      included_features=[col for col in full_dataset_df.columns if col not in del_y] 
      included_features.append(nom_col)
      full_dataset_df=full_dataset_df[included_features] 

   X_train, X_test, y_train, y_test=modelmaker.split_full_dataset_df(full_dataset_df,nom_col, del_y)
   if X_train.empty or X_test.empty or y_train.empty or y_test.empty:
      return None
   return X_train, X_test, y_train, y_test, nom_col 


def models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name, run_dir=None):
   try:
      visualizer.ROCAUC_curve(X_train, y_train, X_test, y_test, RMSD_mode, model,"ROCAUC.png", run_name=run_name,run_dir=run_dir)
   except:
      pass
   try:
      visualizer.precision_recall_curve(X_train, y_train, X_test, y_test, RMSD_mode, model,"precision_recall.png", run_name=run_name,run_dir=run_dir)
   except:
      print("precision_recall failed")
   visualizer.learning_curve(X_train, y_train, model,"learning_curve.png", run_name=run_name,run_dir=run_dir)
   visualizer.feature_importance(X_train, y_train,model,"feature_importance.png", run_name=run_name,run_dir=run_dir)
   visualizer.recursive_feature_elimination(X_train, y_train,model,"recursive_feature_elimination.png", run_name=run_name,run_dir=run_dir)


def dataset_preprocess_viz(full_dataset_df, RMSD_mode, bin_thresholds, val_dic,RMSD_name,run_dir):
   #Visualize Class Distributions
   visualizer.class_distribution(full_dataset_df, RMSD_mode, bin_thresholds, RMSD_name+"_class_distributions.png",  run_dir=run_dir)
   visualizer.num_distribution(full_dataset_df, RMSD_mode,bin_thresholds,val_dic, RMSD_name,RMSD_name+"_distribution.png", run_dir=run_dir)

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
def baseline_model_experimenter(full_dataset_df, RMSD_mode, del_y, models, models_torun_results, full_dataset_df_sele, save_prefix,run_dir=None):
   #TEST PERFORMANCE OF BASELINE MODELS, RUDIMENTARY FEATURES
   included_features=["length", "identity", "similarity"] 
   X_train, X_test, y_train, y_test, nom_col=prepare_for_train(full_dataset_df, RMSD_mode, del_y, included_features=included_features)
   #no hyperparameter optimization
   models_torun_baseline_rud_row, model=baseline_models(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode,run_name="baseline_rud",  run_dir=run_dir)
   models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
   print(models_torun_results)
   models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_rud_non_opti"+save_prefix+".csv"))
   models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name="baseline_rud" ,run_dir=run_dir)

   #with hyperparameter optimization
   models_torun_baseline_rud_row, model=baseline_models_opti(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode,run_name="baseline_rud_opti",run_dir=run_dir)
   models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
   print(models_torun_results)
   models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_rud_opti"+save_prefix+".csv"))
   models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name="baseline_rud_opti",run_dir=run_dir)

   #TEST PERFORMANCE OF BASELINE MODELS, ALL FEATURES
   X_train, X_test, y_train, y_test, nom_col=prepare_for_train(full_dataset_df_sele, RMSD_mode, del_y)
   #no hyperparameter optimization
   models_torun_baseline_rud_row, model=baseline_models(X_train, X_test, y_train, y_test, nom_col, models,RMSD_mode,run_name="baseline_all",run_dir=run_dir)
   models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
   print(models_torun_results)
   models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_all_non_opti"+save_prefix+".csv"))
   models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name="baseline_all",run_dir=run_dir)

   #with hyperparameter optimization
   models_torun_baseline_rud_row, model=baseline_models_opti(X_train, X_test, y_train, y_test, nom_col, models, RMSD_mode, run_name="baseline_all_opti",run_dir=run_dir)
   models_torun_results=pd.concat([models_torun_results, models_torun_baseline_rud_row])
   print(models_torun_results)
   models_torun_results.to_csv(os.path.join(run_dir,"RESULTS_all_opti"+save_prefix+".csv"))
   models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, run_name="baseline_all_opti",run_dir=run_dir)

   return models_torun_results


def full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered, 
                           nr_of_balanced_bins, balanced=False, full_dataset_path=None,fast_mode=None, staggered=False, 
                           stagger_params_all=None, manual_sele=None, save_prefix=None,run_dir=None):
   
   #prepare full dataset with chosen target-y and all features. Also extract validation set, to be used for final testing. 
   with open(os.path.join(run_dir,"info.txt"), 'a') as file:
      file.write("""START OF RUN\nthreshold_bins_staggered: {},\nstep_size_nom_staggered: {},\nnr_of_balanced_bins: {},\nbalanced: {},\nfull_dataset_path: {},\nfast_mode: {},\nsave_prefix: {},\nrun_dir: {}\n""".format(threshold_bins_staggered,step_size_nom_staggered, 
                           nr_of_balanced_bins, balanced, full_dataset_path,fast_mode, save_prefix,run_dir))

   validation_dataset_df, full_dataset_df, complete_df, bin_thresholds, val_dic=prepare_full_dataset_df(model_dir, log_dir, actual_dir, threshold_bins_staggered,
                                                                                 step_size_nom_staggered, nr_of_balanced_bins, full_dataset_path=full_dataset_path, 
                                                                                 fast_mode=fast_mode, save_prefix=save_prefix, run_dir=run_dir, balanced=balanced, 
                                                                                 staggered=staggered, stagger_params_all=stagger_params_all)
   
   #define target_ys
   del_y, local_CA_y_col_list, local_AA_y_col_list,global_CA_y_col_list, global_AA_y_col_list,num_y_col=modelmaker.define_RMSD_modes(full_dataset_df)

   #define the model types used 
   models=modelmaker.all_classifier_types(model_list=None)
   
   #start of testing suite:
   for RMSD_mode, RMSD_name, threshold_bins_staggered_RMSD, step_size_nom_staggered_RMSD in zip([local_CA_y_col_list, local_AA_y_col_list,global_CA_y_col_list, global_AA_y_col_list], num_y_col, threshold_bins_staggered, step_size_nom_staggered):
      if RMSD_mode!=global_CA_y_col_list:
         continue
      #feature selection
      dataset_preprocess_viz(complete_df, RMSD_mode, bin_thresholds, val_dic,RMSD_name,run_dir=run_dir)
      with open(os.path.join(run_dir,"info.txt"), 'a') as file:
         file.write("FEATURE SELECTION\n")
         file.write("Feature selection with max number of features : {}\n".format(30))
      selection_dic=preprocessor.feature_selection(full_dataset_df, RMSD_mode, 30, run_dir, manual_sele=manual_sele)
      run_dir_og=run_dir
      full_dataset_df_og=full_dataset_df
      for selection in selection_dic.keys():
         run_dir=run_dir_og
         full_dataset_df=full_dataset_df_og
         all_cols=np.concatenate([selection_dic[selection],del_y],axis=0)
         full_dataset_df_sele=full_dataset_df[all_cols]
         run_dir=os.path.join(run_dir,str(selection))
         if not os.path.isdir(run_dir):
            os.makedirs(run_dir)
         with open(os.path.join(run_dir,"info.txt"), 'a') as file:
            file.write("STARTING RUN WITH SELECTION\n".format(selection))
         with open(os.path.join(run_dir,"selected_features.txt"), 'w') as file:
            file.writelines(selection_dic[selection])
         vif=preprocessor.calc_vif_score(full_dataset_df[selection_dic[selection]], del_y)
         vif.to_csv(os.path.join(run_dir,"vif_features_common_sele.csv"))
         with open(os.path.join(run_dir,"info.txt"), 'a') as file:
            file.write("STARTING MODEL GENERATION FOR RMSD MODE: {}\n".format(RMSD_mode))

         #make dataset with only the current RMSD_mode's target y columns
         RMSD_cols=[col for col in full_dataset_df_sele if col in del_y and col not in RMSD_mode]
         full_dataset_df_RMSD_sele=full_dataset_df_sele[full_dataset_df_sele.columns[~full_dataset_df_sele.columns.isin(RMSD_cols)]] #remove all target cols not relevant for this RMSD mode
         full_dataset_df=modelmaker.clean_dataset(full_dataset_df_RMSD_sele)
        
         #MODEL MAKER FULL RUN#######################################################################################################
         models_torun_results=pd.DataFrame() #dataframe containing final MCC of the best model for each specified run
         for column in full_dataset_df_RMSD_sele.columns: ################TEMPORARY
            try:
               full_dataset_df_RMSD_sele[column]=full_dataset_df_RMSD_sele[column].to_numpy().astype(np.float)
            except:
               del full_dataset_df_RMSD_sele[column]
               with open(os.path.join(run_dir,"info.txt"), 'a') as file:
                  file.write("Deleted the following column as it was not encoded: {}\n".format(column))
         ################
         #BASELINE MODELS
         ################
         #models_torun_results=baseline_model_experimenter(full_dataset_df_RMSD_sele, RMSD_mode, del_y, models, models_torun_results, 
         #                     full_dataset_df_RMSD_sele, save_prefix,  run_dir=run_dir)
         
         ################
         #STACKED MODELS
         ################
         with open(os.path.join(run_dir,"info.txt"), 'a') as file:
            file.write("STARTING STACKED MODEL\n")
         for target_y in RMSD_mode:
            if "nom" in target_y:
               print("nominal classes: ",np.unique(full_dataset_df_RMSD_sele[target_y]))
         models_torun_double_staggered, model, X_train, X_test, y_train, y_test, MCC_of_bins=modelmaker.double_staggered_run(full_dataset_df_RMSD_sele, del_y,models, RMSD_mode, val_dic, run_dir=run_dir)
         models_torun_results=pd.concat([models_torun_results, models_torun_double_staggered, MCC_of_bins])
         models_torun_results.to_csv(os.path.join(run_dir,"OVERNIGHT_RESULTS_stacked"+save_prefix+".csv"))

         models_viz(X_train, y_train, X_test, y_test, RMSD_mode, model, "stacked_models"+save_prefix, run_dir=run_dir)




if __name__ == '__main__':
   model_dir=sys.argv[1]
   log_dir=sys.argv[2]
   actual_dir=sys.argv[3]
   threshold_bins_staggered=[[1,7,8], [1,8,6],[2,10,8],[1,10,8]] #"local_AA", "local_CA", "global_AA", "global_CA"
   step_size_nom_staggered=[2,2,2,2]
   nr_of_balanced_bins=10

   manual_sele=[['tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff'], 
   ['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', "contacts_all"], 
   ['blosum_dist', 'contacts_all', 'Screlacc', 'Access', 'simlength', 'Scacc',  'identity', 'Relacc', 'tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff'], 
   ['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', 'simlength', "PCAtip_res_blosum_621_blosum_62_","PCAtip_res_blosum_622_blosum_62_","PCAtip_res_blosum_623_blosum_62_"],
   ['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', 'simlength', "PCAtip_res_NLF1_NLF_col","PCAtip_res_NLF2_NLF_col","PCAtip_res_NLF3_NLF_col"],
   ['blosum_dist','tip_pos', 'protrusion', 'length', 'total_charge', 'nr_charged', 'identity', 'similarity','Hydropathy','Hydropathy_diff', 'simlength', "res_charge_physiochem_4d_","res_sc_nr_physiochem_4d_","res_compactness_physiochem_4d_"]]


   """full_dataset_path="/serv/www/html_lilian/CDRH3loop/CDRH3loop/staggered_try_out/staggered_try_out.csv"
   bin_thresholds=[1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9]
   val_dic={"nom_label":[1.0,2.0,3.0], "lower_threshold":[0,2,4], "upper_threshold":[2,4,100]}
   full_dataset_pack=[full_dataset_path, bin_thresholds, val_dic]
   save_prefix="staggered_0.2_5_1.5_2"
   stagger_params_all=[[0.2,6,1.5,2],[0.2,5,1.5,2],[0.2,9,1.5,2],[0.2,9,1.5,2]] #[staggering_value, max_angstrom, first_layer_size, second_layer_size]
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, staggered=True,stagger_params_all=stagger_params_all,
                        full_dataset_path=full_dataset_pack, manual_sele=manual_sele,save_prefix=save_prefix, run_dir=run_dir)
"""
   full_dataset_path="/serv/www/html_lilian/CDRH3loop/CDRH3loop/staggered_0.1_5_1_2_global_CA/staggered_0.1_5_1_2_global_CA.csv"
   bin_thresholds=[1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,
   3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,
   7.5,7.6,7.7,7.8,7.9,8.0,8.1,8.2,8.3,8.4,8.5,8.6,8.7,8.8,8.9]
   val_dic={"nom_label":[1.0,2.0,3.0], "lower_threshold":[0,2,4,6,8], "upper_threshold":[2,4,6,8,100]}
   full_dataset_pack=[full_dataset_path, bin_thresholds, val_dic]
   
   save_prefix="staggered_0.1_5_1_2_global_CA"
   stagger_params_all=[[0.1,6,1.5,2],[0.1,5,1,2],[0.1,9,1,2],[0.1,9,1,2]] #[staggering_value, max_angstrom, first_layer_size, second_layer_size]
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, staggered=True,stagger_params_all=stagger_params_all,
                        full_dataset_path=full_dataset_pack, manual_sele=manual_sele,save_prefix=save_prefix, run_dir=run_dir)
                        
   full_dataset_pack=None
   save_prefix="staggered_0.1_5_1.5_1"
   stagger_params_all=[[0.1,6,1.5,1],[0.1,5,1.5,1],[0.1,9,1.5,1],[0.1,9,1.5,1]] #[staggering_value, max_angstrom, first_layer_size, second_layer_size]
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, staggered=True,stagger_params_all=stagger_params_all,
                        full_dataset_path=None, manual_sele=manual_sele,save_prefix=save_prefix, run_dir=run_dir)
   

   full_dataset_path=None
   threshold_bins_staggered=[[1,7,8], [1,8,6],[2,10,8],[1,10,8]] #"local_AA", "local_CA", "global_AA", "global_CA"
   step_size_nom_staggered=[2,2,2,2]
   nr_of_balanced_bins=10
   save_prefix="186_2_uniform"
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, balanced=False,
                        full_dataset_path=full_dataset_path, manual_sele=manual_sele, save_prefix=save_prefix, run_dir=run_dir)

   full_dataset_path=None
   threshold_bins_staggered=[[1,7,9], [1,8,12],[2,10,9],[1,10,9]] #"local_AA", "local_CA", "global_AA", "global_CA"
   step_size_nom_staggered=[3,3,3,3]
   nr_of_balanced_bins=10
   save_prefix="1812_3_uniform"
 
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, balanced=False,
                        full_dataset_path=full_dataset_path, save_prefix=save_prefix, run_dir=run_dir)

   full_dataset_path=None
   threshold_bins_staggered=[[1,7,9], [1,8,12],[2,10,9],[1,10,9]] #"local_AA", "local_CA", "global_AA", "global_CA"
   step_size_nom_staggered=[4,4,4,4]
   nr_of_balanced_bins=10
   save_prefix="1812_4_uniform"
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, balanced=False,
                        full_dataset_path=full_dataset_path, save_prefix=save_prefix, run_dir=run_dir)

   full_dataset_path=None
   threshold_bins_staggered=[[1,7,8], [1,8,6],[2,10,8],[1,10,8]] #"local_AA", "local_CA", "global_AA", "global_CA"
   step_size_nom_staggered=[2,2,2,2]
   nr_of_balanced_bins=10
   save_prefix="186_2_balanced"
   #now = datetime.now()
   #save_prefix=now.strftime("%d%b_%H%M%S")+save_prefix
   run_dir=os.path.join('/serv/www/html_lilian/CDRH3loop/CDRH3loop', save_prefix)
   if not os.path.isdir(run_dir):
      os.makedirs(run_dir)
   full_experimentation(model_dir, log_dir, actual_dir, threshold_bins_staggered,step_size_nom_staggered,nr_of_balanced_bins, balanced=True,
                        full_dataset_path=full_dataset_path, save_prefix=save_prefix, run_dir=run_dir)
