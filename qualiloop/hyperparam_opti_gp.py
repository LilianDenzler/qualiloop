#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Hyperparameter Optimizer- Gaussian Minimizer
   File:       modelmaker.py
   
   Version:    V1.1
   Date:       14.03.21
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
   The program contains all functions needed for hyperparameter optimization
   of the models used for evaluation or for building the final model. 
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
import inspect
import numpy as np


from skopt.space import Real, Integer, Categorical
from skopt.utils import use_named_args
from sklearn.metrics import make_scorer

from sklearn.model_selection import cross_val_score
from sklearn.metrics import log_loss
from scipy.optimize import minimize
from sklearn.metrics import matthews_corrcoef
from sklearn.exceptions import ConvergenceWarning, FitFailedWarning
from sklearn.utils._testing import ignore_warnings 


def toplayer(models,X_train, y_train, X_test, y_test):
   model_dic={} #add all models here
   model_list=[]
   for model_name,model in models:
      model_list.append(model)

   for model_name, model in models:
     
      index=models.index((model_name, model))
      print(model_name)
      model_dic=param_dics(model_name,model_list,models)
      space=model_dic
      optimize_gp(space)
      print(grid_best_params)
      model=model.set_params(**grid_best_params)
      models[index]=(model_name,model)

      
      if model_name=="Voting_all_hard" or model_name=="Voting_all_soft":
         param_grid=define_grid_vote(model_dic,models,X_test, X_train, y_train, y_test)
         model=model.set_params(**param_grid)
         models[index]=(model_name,model)
      
   return models

######################################
#!!!!!!!!!!!!!define ranges in dics, dont use some int etc
def param_dics(model_name, model_list, models):
   #define all the parameter dictionary grids
   ################################################################################################################################


   LogisticRegression_dic=[
      Categorical(['l1', 'l2', 'elasticnet', 'none'],name='penalty'),
      Real(100,100**-4, "log-uniform", name='C'),
      Categorical(['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'],name='solver'),
   ]

   LinearDiscriminantAnalysis_dic={
      'solver': ['svd', 'lsqr', 'eigen'],
      'shrinkage':['auto',dists.uniform(0, 1)],
      #'tol': tuple(np.linspace(0, 500, 100))
   }

   KNeighborsClassifier_dic={
      'n_neighbors':dists.randint(0,500),
      'weights': ['uniform', 'distance'],
      'algorithm':['auto', 'ball_tree', 'kd_tree', 'brute'],
      'leaf_size': dists.randint(0,500),
   }

   DecisionTreeClassifier_dic={
      'criterion': ['gini', 'entropy'],
      'splitter': ['best', 'random'],
      #'max_depth': [dists.randint(0,500),None],
      'min_samples_split': dists.uniform(0, 1), 
      #'min_samples_leaf': dists.uniform(0, 0.5),
      'min_weight_fraction_leaf': dists.uniform(0, 0.5),
      #'max_features': ['auto', 'sqrt', 'log2',dists.uniform(0, 1)],
      'random_state':[None, 42],
      #'max_leaf_nodes': [dists.randint(0,500),None],
      #'min_impurity_decrease': dists.uniform(0, 500),
      #'min_impurity_split': dists.uniform(0, 500),
      'class_weight': ["balanced", None],
      #'ccp_alpha': dists.uniform(0, 500),
   }

   GaussianNB_dic={
      #'var_smoothing': dists.uniform(0, 0.5)
   }

   RandomForestClassifier_dic={
      'n_estimators':dists.randint(0,500),
      'criterion': ['gini', 'entropy'],
      #'max_depth': [dists.randint(0,500),None],
      'min_samples_split':dists.uniform(0, 1),
      #'min_samples_leaf': dists.uniform(0, 0.5),
      #'min_weight_fraction_leaf': dists.randint(0,0.5),
      #'max_features':['auto', 'sqrt', 'log2',dists.uniform(0, 1)],
      #'max_leaf_nodes': [dists.randint(0,500),None],
      #'min_impurity_decrease': dists.uniform(0, 500),
      #'min_impurity_split': dists.uniform(0, 500),
      'bootstrap': [True, False],
      'oob_score': [True, False],
      'random_state': [None, 42],
      'class_weight': ['balanced', 'balanced_subsample', None],
      #'ccp_alpha': dists.uniform(0, 500),
      #'max_samples': [dists.uniform(0, 1),None]
   }

   SVC_dic={
      'C': dists.uniform(0, 500),
      'kernel': ['linear', 'poly', 'rbf', 'sigmoid'],
      'degree': dists.randint(0,500),
      'gamma': ['scale', 'auto',dists.uniform(0, 500)],
      'coef0': dists.uniform(0, 500),
      'shrinking': [True, False],
      'probability': [True, True],
      #'tol': tuple(np.linspace(0, 500, 100)),
      'cache_size': dists.uniform(0, 500),
      'class_weight': ['balanced', None],
      'max_iter': dists.randint(0,500),
      'decision_function_shape': ['ovo', 'ovr'],
      #'break_ties': [True,False],
      'random_state': [None, 42]
   }

   Voting_all_hard_dic={
      'estimators':[models,models],
      'weights': [[1]*len(model_list), None]
   }

   Voting_all_soft_dic={
      'estimators':[models,models],
      'weights': [[1]*len(model_list), None]
   }

   association_dic={
      "LogisticRegression":LogisticRegression_dic, 
      "LinearDiscriminantAnalysis": LinearDiscriminantAnalysis_dic, 
      "KNeighborsClassifier": KNeighborsClassifier_dic, 
      "DecisionTreeClassifier": DecisionTreeClassifier_dic, 
      "GaussianNB": GaussianNB_dic, 
      "RandomForestClassifier": RandomForestClassifier_dic,
      "SVC": SVC_dic,
      "Voting_all_soft":Voting_all_soft_dic,
      "Voting_all_hard":Voting_all_hard_dic
   }
   model_dic=association_dic[model_name]

   return model_dic
@use_named_args(space)
def loss_func(**params) :
   """
   Determine the optimal hyperparameters of the model using a random seaching method. 


   Input:  model --- the model for which the hyperparameters are being optimized
         X_train ---  X feature training set
         y_train  ---  training set of the target y value 
         X_test  ---   X feature testing set for hyperparameter optimizatio
         y_test  ---  testing set of the target y value for hyperparameter optimization

    Return: dictionary of all optimal hyperparameters for the model
         e.g. {'bootstrap': True,
             'max_depth': 70,
             'max_features': 'auto',
             'min_samples_leaf': 4,
             'min_samples_split': 10,
             'n_estimators': 400}
    """
   #make MCC a scorable metric -> the hyperparameters will be chose so that MCC  is maximized, not the accuracy as is default 
   model.set_params(**params)
   mcc_scorer=make_scorer(matthews_corrcoef)
   model.fit(X_train, y_train)
   res=np.mean(cross_val_score(model,X_train, y_train, cv=10, n_jobs=-1, scoring=mcc_scorer))
   print(model, res)
   return res

def optimize_gp(space):
   model_gp=gp_minimize(loss_func,space, n_calls=50, random_state=0)
   print(model_gp.fun)
   print(model_gp.x)
   quit()


def log_loss_func(weights,predictions, y_test):
   ''' scipy minimize will pass the weights as a numpy array '''
   final_prediction = 0
   for weight, prediction in zip(weights, predictions):
      #print(prediction)
      final_prediction += weight*prediction

   return log_loss(y_test, final_prediction)
       
      
def optimize_weights(models,X_test, X_train, y_train, y_test):
   #https://www.kaggle.com/hsperr/finding-ensamble-weights
   #print("MODELS",models,"MODELS")
   predictions = []
   for model_name, model in models:
      #print(model_name)
      model.fit(X_train, y_train)
      model_pred=model.predict_proba(X_test)
      #print(model_name,model_pred)
      predictions.append(model_pred)
   starting_values = [0.5]*len(predictions)

   #adding constraints  and a different solver as suggested by user 16universe
   #https://kaggle2.blob.core.windows.net/forum-message-attachments/75655/2393/otto%20model%20weights.pdf?sv=2012-02-12&se=2015-05-03T21%3A22%3A17Z&sr=b&sp=r&sig=rkeA7EJC%2BiQ%2FJ%2BcMpcA4lYQLFh6ubNqs2XAkGtFsAv0%3D
   cons = ({'type':'eq','fun':lambda w: 1-sum(w)})
   #our weights are bound between 0 and 1
   bounds = [(0,1)]*len(predictions)
   #print(predictions)
   res = minimize(log_loss_func,starting_values, args=(predictions,y_test),method='SLSQP', bounds=bounds, constraints=cons) # also try method='Nelder-Mead'
   weights=res['x']

   return weights

def define_grid_vote(model_dic,models,X_test, X_train, y_train, y_test):
   print("define_grid_vote")
   tup_dict = dict(models) 
   tup_dict.pop('Voting_all_hard')
   tup_dict.pop('Voting_all_soft')
   models = tuple(tup_dict.items())
   for key in model_dic.keys():
      print("KEY", key)
      if key=='weights':
         print("INLOOP")
         param_space=optimize_weights(models,X_test, X_train, y_train, y_test)
         print("PARAM_SPACE",param_space)
      
   weights= param_space
   #print(model_dic)
   model_dic['weights']=weights
   model_dic['estimators']=models

   return model_dic

