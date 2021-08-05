#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Final Model
   File:       final_predictor.py
   
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
   The program contains the fucntions needed to get a prediction form the 
   finalized model for a previously unseen instance.
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
from qualiloop import myfunctions
from qualiloop import visualizer
import model_experimenter

import pickle
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import VarianceThreshold
import xgboost as xgb
import tensorflow as tf

import doctest
from qualiloop import preprocessor

import argparse
import configparser
from datetime import datetime
from pathlib import Path
import glob
import keras
from keras.models import load_model
import ast

from bson.binary import Binary
from bson import ObjectId
from pymongo import MongoClient
os.environ["CONFIG"] = "/home/lilian/sync_project/WWW/CDRH3loop/CDRH3loop/qualiloop/final_config.conf"


def get_file_paths(touse_model_dir):
   val_df_path=os.path.join(touse_model_dir,'val_df.csv')
   if os.path.exists(val_df_path)==True:
      print(val_df_path)
   else:
      val_df_path=glob.glob(os.path.join(touse_model_dir,'/**/val_df.csv'), recursive=True)
      val_df_path=val_df_path[0]
      #print('**/matrices_blosum_pam/'+matrix)
      if len(val_df_path)>1:
         raise Warning('more than one val_df file was found')
      elif len(val_df_path)==0:
         raise Warning("val_df_path file is not found!")

   selected_features_path=os.path.join(touse_model_dir,'selected_features.csv')
   if os.path.exists(selected_features_path)==True:
      print(selected_features_path)
   else:
      selected_features_path=glob.glob(os.path.join(touse_model_dir,'/**/selected_features.csv'), recursive=True)
      selected_features_path=selected_features_path[0]
      if len(selected_features_path)>1:
         raise Warning('more than one selected_features file was found')
      elif len(val_df_path)==0:
         raise Warning("selected_features file is not found!")

   return val_df_path, selected_features_path

def get_val_dic(touse_model_dir):
   val_df_path,selected_features_path=get_file_paths(touse_model_dir)
   print(val_df_path)
   val_dic=pd.read_csv(val_df_path,index_col=0, header=None)
   print(val_dic)
   threshlow=val_dic.loc["lower_threshold",1]
   print(threshlow, type(threshlow))
   threshlow =np.array(threshlow[1:-1].split(),float)
   print(type(threshlow))
   print(threshlow)
   
   threshup=val_dic.loc["upper_threshold",1]
   threshup =np.array(threshup[1:-1].split(),float)
   bin_thresholds=val_dic.loc["bin_thresholds",1]
   bin_thresholds =np.array(bin_thresholds[1:-1].split(),float)
   

   #get selected features
   selection=pd.read_csv(selected_features_path, index_col=None, header=0)
   selection=selection[["sele"]].to_numpy().tolist()
   selection= [item for sublist in selection for item in sublist]
   print(selection)

   return threshlow,threshup,bin_thresholds,selection

def make_full_df(model_dir, log_dir):
   full_dataset_df=pd.DataFrame()
   for file in os.listdir(model_dir):
      print(model_dir, file)
      file=os.path.join(model_dir, file)
      filename=os.path.splitext(os.path.basename(file))[0]
      filename=filename.replace(".pdb", "")
      log_name=os.path.join(log_dir, filename+".log")
      if os.path.isfile(file)==False or os.path.isfile(log_name)==False:
         print(f"either {file}, or {log_name} could not be found. The model named {filename} will be excluded. ")
         continue
      elif os.path.getsize(file)==0 or os.path.getsize(log_name)==0:
         print(f"either {file}, or {log_name} are empty files. The model named {filename} will be excluded. ")
         continue
      all_features_df=modelmaker.more_features(file,log_name)
      all_features_df.reset_index(drop=True, inplace=True)
      full_dataset_df.reset_index(drop=True, inplace=True)
      full_dataset_df=pd.concat([full_dataset_df, all_features_df])
   return full_dataset_df


def retrain_prep(model_dir, actual_dir, log_dir, retrain_df_path=None):
   if retrain_df_path:
      retrain_df=pd.read_csv(os.path.join(retrain_df_path),header=0)
   else:
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

            RMSD_num_df=modelmaker.get_RMSD_num(file, actual_file, fast_mode=None)
            ID_df=pd.DataFrame(np.array([filename]), columns=["ID"])
            all_features_df=modelmaker.more_features(file,log_name)

            answers_df=pd.concat([ID_df,RMSD_num_df],axis=1)
            all_features_df=pd.concat([all_features_df,ID_df], axis=1)
            full_dataset_df=pd.concat([all_features_df,full_dataset_df],axis=0)
            full_dataset_df.reset_index(drop=True, inplace=True)
            full_answers=pd.concat([full_answers, answers_df],axis=0)
            full_answers.reset_index(drop=True, inplace=True)
            answer_col=full_answers
   return full_dataset_df


def use_models(RMSD_mode, threshup, threshlow, selection, touse_model_dir, full_dataset_df, mode):
   print(full_dataset_df[selection].columns)
   print(len(full_dataset_df[selection].columns))

   for target_y in RMSD_mode:
      try:
         print(full_dataset_df.columns)
         best_model_of_bin=xgb.XGBClassifier()
         best_model_of_bin.load_model(os.path.join(touse_model_dir, "layer_1","bin"+str(target_y)+".json"))
         print(selection)
         print(full_dataset_df.columns)
         sele_df=full_dataset_df[selection]
         print(sele_df)
         prediction = best_model_of_bin.predict_proba(full_dataset_df[selection])
         print(prediction)
         
      except:
         best_model_of_bin=joblib.load(os.path.join(touse_model_dir, "layer_1","bin"+str(target_y)+".pkl"))
         sele_df=full_dataset_df[selection]
         print(sele_df)
         prediction = best_model_of_bin.predict_proba(full_dataset_df[selection])
         print(prediction)
      full_dataset_df[target_y]=prediction[:,1]

   predicted_cols = [col for col in full_dataset_df.columns if "bin" in col]
   predicted_df=full_dataset_df[predicted_cols]
   second_layer_df=pd.DataFrame()
   for i, i_low in zip(threshup, threshlow):
      if i==100:
         continue
      multiply_list=[]
      for key,value in predicted_df.items():
         key=key.split("bin", 1)[1]
         key=float(key)
         if key <=i and key>i_low:
            multiply_list.append([value])
      number_of_vals=len(multiply_list)
      arr = np.array(multiply_list)
      sum_array=arr.sum(axis=0)/number_of_vals
      print("SUMM ARRAY", sum_array)

      sum_array=sum_array[0]   
      print("SUMM ARRAY", sum_array)
      second_layer_df[str("nom")+str(i)]=sum_array
   
      print("test_df",second_layer_df)
      #second_layer_df.columns = second_layer_df.columns.astype(str)
      full_dataset_df=full_dataset_df[[col for col in full_dataset_df.columns if col not in predicted_cols]]
      target_y=mode+"_nom"
      #full_dataset_df=pd.concat([full_dataset_df, second_layer_df],axis=1)
      
   for model_mode in ("all","only_sums"):
      if mode=="all":
         input(selection)
         sele_df=full_dataset_df[selection]
         print(sele_df)
         X=pd.concat([sele_df, second_layer_df],axis=1)
      else:
         X=second_layer_df
      if os.path.exists(os.path.join(touse_model_dir, model_mode +".json")) or os.path.exists(os.path.join(touse_model_dir, model_mode +".pkl")):
         try:
            print(sele_df.columns)
            best_model_second_layer=xgb.XGBClassifier()
            best_model_second_layer.load_model(os.path.join(touse_model_dir, model_mode +".json"))
            prediction = best_model_second_layer.predict(X)
            full_dataset_df["prediction_"+model_mode]=prediction
            print(prediction)
            print(full_dataset_df)
            answer_df=pd.concat([full_dataset_df],axis=1)
               
         except:
            best_model_second_layer=joblib.load(os.path.join(touse_model_dir, model_mode +".pkl"))
            prediction = best_model_second_layer.predict(X)
            print(prediction.shape,full_dataset_df.shape)
            full_dataset_df["prediction_"+model_mode]=prediction
            print(prediction)
            print(full_dataset_df)
            answer_df=pd.concat([full_dataset_df],axis=1)
      else:
         print(f"NO MODEL WAS FOUND AT {os.path.join(touse_model_dir,'bin_models',model_mode)}")
   print(answer_df)
   return answer_df

def tilted_loss(q,y,f):
   e = (y-f)
   return K.mean(K.maximum(q*e, (q-1)*e), axis=-1)

def use_numerical(touse_model_dir,X_test):
   try:
      print(X_test)
      model=joblib.load(os.path.join(touse_model_dir))
      y_pred = model.predict(X_test)
   except:
      print(touse_model_dir)
      model = keras.models.load_model(touse_model_dir)
      X_test_tf = tf.convert_to_tensor(X_test)
      try:
         model = Sequential()
         model.add(Dense(units=3, input_dim=X_train.shape[1],activation='relu'))
         model.add(Dense(units=3, activation='relu'))
         model.add(Dense(1))
         model.compile(loss='mae', optimizer='adadelta')
         model.load_weights(touse_model_dir)
         y_pred = model.predict(X_test_tf)
      except:
        return("Could not load model correctly.")
      y_pred = model.predict(X_test_tf)
      y_pred= y_pred.flatten()

   """model = XGBRegressor(n_estimators=1000, max_depth=7, eta=0.1, subsample=0.7, colsample_bytree=0.8, reg_lambda=0.8, gamma=24, learning_rate=0.123)
   if already_trained==True:
      model.load_model(os.path.join("XGBRegressor"+".json"))"""
   input(y_pred)
   return y_pred


def run_with_unknown(mode,model_type,touse_model_dir, ml_model_path, model_dir, log_dir, actual_dir=None, re_train=False, run_dir=None, full_dataset_df=None):
   threshlow,threshup, bin_thresholds,selection=get_val_dic(touse_model_dir)
   RMSD_mode=[mode+"_bin"+str(i) for i in bin_thresholds]
   print(RMSD_mode)
   print(model_type)
   #if model is not to be re-trained, the model dir and log dir files are used to create a dataframe with all features
   #if model is being re-trained, the dataframe will also contain the RMSD values, as well as bin and nom RMSD values in addition to features
   if full_dataset_df:
      full_dataset_df=pd.read_csv(os.path.join(full_dataset_df),header=0, index_col=False)
   elif re_train==False:
      full_dataset_df=make_full_df(model_dir, log_dir)
   elif re_train==True:
      full_dataset_df=retrain_prep(model_dir, actual_dir, log_dir)
   full_dataset_df.to_csv(os.path.join(run_dir,"full_dataset_df.csv"), index=None)
   full_dataset_df=preprocessor.scaling(full_dataset_df,"robust", run_dir)
   full_dataset_df=model_experimenter.encode_full_df(full_dataset_df, run_dir)
   full_dataset_df.to_csv(os.path.join(run_dir,"full_dataset_df_encoded.csv"), index=None)
   print(full_dataset_df.columns)
   full_dataset_df=full_dataset_df[selection]
   print(full_dataset_df.columns)

   if model_type=="staggered" or model_type=="equal" or model_type=="balanced":
      if re_train==False:
         prediction_df=use_models(RMSD_mode, threshup, threshlow, selection, ml_model_path, full_dataset_df,mode)
      elif re_train==True:
         retrain_models(RMSD_mode, threshup, threshlow, selection, ml_model_path, full_dataset_df)

   if model_type=="numerical":
      y_pred=use_numerical(ml_model_path,full_dataset_df)

   quit()


def param_parser():
   """ {'name': 'hiqw', 'model_type': 'equal', 'mode': 'local_CA', 
   'sele': 'manual', 'feature': ['loop_seq', 'loop_length', 
   'total_charge', 'nr_charged', 'tip_res', 'protrusion', 'nr_sad', 
   'happiness_mean', 'template']}"""

   conf_parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,add_help=False)
   conf_parser.add_argument("-config",dest="config", help="Specify config file", metavar="FILE")
   args, remaining_argv = conf_parser.parse_known_args()
   print("args",args)

   settings = {}

   if args.config:
      config = configparser.SafeConfigParser()
      config.read([args.config])
      print(config.items("model_1"))
      settings.update(dict(config.items("model_1")))
   input(settings)

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

   saved_classifier = parser.add_mutually_exclusive_group()
   saved_classifier.add_argument('-classifier' , metavar='classifier_keyword',choices=['model_1', 'model_2', 'model_3'], help='',default=None)
   saved_classifier.add_argument('-cp', nargs='?',type=Path, metavar="classifierdirectory ", help='path of directory where all classifier files are stored')

   parser.add_argument('-name',default=None, type=Path,help='save the project run under this name.')
   parser.add_argument('-r',  default=None,type=Path,help='save your project in this directory.')
   parser.add_argument('-type' , metavar='type',choices=['staggered', 'equal', 'balanced', 'numerical'], help='',default=None)
   parser.add_argument('-mode',type=str,metavar="mode", help='', default=None)
   parser.add_argument('-ml', nargs='?',type=Path, metavar="classifier path", help='path of actual ML model')

   parser.add_argument('-model_dir',type=Path, metavar="model dir path", required=True, help='path of actual ML model')
   parser.add_argument('-log_dir',type=Path, metavar="log dir path", required=True,help='path of actual ML model')

   parser.add_argument('-full',type=Path, metavar="full dataset df",default=None,help='path of actual ML model')


   
   #parser.add_argument('-config',metavar="FILE", dest="config_path", help='path to config file.')
   args2 = vars(parser.parse_args())
   print(args)
   input(args2)

   settings.update({k: v for k, v in args2.items() if v is not None})
   settings.update({k: v for k, v in args2.items() if v is None and k not in settings.keys()}) 
   for key,val in settings.items():
      try:
         val=eval(val)
         settings.update({key:val})
      except:
         pass
   if settings["classifier"] ==None and settings["cp"]==None:
      parser.error("Must specify model using either -classifier or -cp")

   if settings["type"]==None:
      parser.error("the following arguments is required: type")
   if settings["mode"]==None:
      parser.error("the following arguments is required: mode")

   
   run_dir="test_final"
   if not os.path.exists(os.path.join(run_dir)):
         os.makedirs(os.path.join(run_dir))

   #open(os.path.join(run_dir,"info.txt"), 'w')

   run_with_unknown(settings["mode"],settings["type"],settings["cp"], settings["ml"],settings["model_dir"], settings["log_dir"], actual_dir=None, re_train=False, run_dir=run_dir, full_dataset_df=settings["full"])

def mongo_retrieve(type, mode, UserID):
   myclient = pymongo.MongoClient("mongodb+srv://m001-student:m001-mongodb-basics@sandbox.quoyx.mongodb.net/test")
   mydb = myclient.models
   data = mycol.find({"type":type, "mode":mode})
   os.makedir(UserID+"mongo_files")
   for i in data:
      with open(i["name"], "wb") as f:
         f.write(i['file'])

def from_website(mode, model_type, model_dir, log_dir, name, sele, features):
   config_file=os.environ["CONFIG"]
   config = configparser.SafeConfigParser()
   config.read([config_file])
   for section in config.sections():
      correct=True
      config_dic=dict(config.items(section))
      print(config_dic["mode"],config_dic["model_type"],config_dic["sele"])
      print("2",mode, model_type, sele)
     
      if config_dic["mode"]==mode and config_dic["model_type"]==model_type and config_dic["sele"]==sele:
         if sele=="manual":
            if config_dic["features"]==features:
               break
         else:
            break
      else:
         correct=False
   if correct==False:
      sys.stderr.write("error, no config file given")
      quit()
   settings={}
   settings.update(config_dic)
   settings.update({"name":name})
   print(settings)
   run_with_unknown(settings["mode"],settings["model_type"],settings["cp"], settings["ml"],model_dir, log_dir, actual_dir=None, re_train=False, run_dir=settings["name"], full_dataset_df=None)
   

if __name__ == '__main__':
   """mode="local_CA"
   touse_model_dir=sys.argv[1]
   model_dir=sys.argv[2]
   log_dir=sys.argv[3]
   actual_dir=sys.argv[4]"""
   param_parser()
