#!/usr/bin/python3
"""/*************************************************************************
   Program:    Qualiloop Hyperparameter Optimizer- Genetic Algorithm
   File:       hyperparam_opti_genetic.py
   
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
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef
import xgboost as xgb
import random
from sklearn.model_selection import train_test_split
import time
from qualiloop import diff_models_train_genetic_opti
from qualiloop import diff_models_initialize_genetic_opti
from sklearn.metrics import log_loss
from scipy.optimize import minimize


def print_defaults(models):
	for model_name, model in models:
		if model_name=="XGBoost":
			continue
		liner = str(model.__doc__).split('Parameters\n    ----------\n')[1].split('\n\n    Attributes\n')[0].replace('\n        ', '\n').splitlines()
		for i in liner:
			if " : " in i: 
				print(i)

def fitness_mccscore(y_true, y_pred):
	fitness=matthews_corrcoef(y_true, y_pred, sample_weight=None ) #maybe round((..),4)
	return fitness#train the data annd find fitness score


#select parents for mating
def new_parents_selection(population, fitness, numParents):
	selectedParents = np.empty((numParents, population.shape[1])) #create an array to store fittest parents
	
	#find the top best performing parents
	for parentId in range(numParents):
		bestFitnessId = np.where(fitness == np.max(fitness))
		bestFitnessId  = bestFitnessId[0][0]
		selectedParents[parentId, :] = population[bestFitnessId, :]
		fitness[bestFitnessId] = -1 #set this value to negative so this parent is not selected again
	return selectedParents

'''
Mate these parents to create children having parameters from these parents (we are using uniform crossover method)
'''
def crossover_uniform(parents, childrenSize):
	
	crossoverPointIndex = np.arange(0, np.uint8(childrenSize[1]), 1, dtype= np.uint8) #get all the index
	crossoverPointIndex1 = np.random.randint(0, np.uint8(childrenSize[1]), np.uint8(childrenSize[1]/2)) # select half  of the indexes randomly
	crossoverPointIndex2 = np.array(list(set(crossoverPointIndex) - set(crossoverPointIndex1))) #select leftover indexes
	
	children = np.empty(childrenSize)
	
	'''
	Create child by choosing parameters from two parents selected using new_parent_selection function. The parameter values
	will be picked from the indexes, which were randomly selected above. 
	'''
	for i in range(childrenSize[0]):
		
		#find parent 1 index 
		parent1_index = i%parents.shape[0]
		#find parent 2 index
		parent2_index = (i+1)%parents.shape[0]
		#insert parameters based on random selected indexes in parent 1
		children[i, crossoverPointIndex1] = parents[parent1_index, crossoverPointIndex1]
		#insert parameters based on random selected indexes in parent 1
		children[i, crossoverPointIndex2] = parents[parent2_index, crossoverPointIndex2]
	return children


###################################################################################
#IMPLEMENT
def implement_opti(X_train,y_train,X_test,y_test, model, model_name, number_of_parameters):
	number_of_parents = 10 #number of parents to start
	number_of_parents_mating = 5 #number of parents that will mate
	number_of_generations = 100 #number of genration that will be created

	populationSize = (number_of_parents, number_of_parameters)
	#initialize the population with randomly generated parameters
	initialize_population_to_call = getattr(diff_models_initialize_genetic_opti, str(model_name)+"_initialize_pop")
	population,param_names, list_encode_dics = initialize_population_to_call(number_of_parents)
	
	#define an array to store the fitness  hitory
	fitnessHistory = np.empty([number_of_generations+1, number_of_parents])
	#define an array to store the value of each parameter for each parent and generation
	populationHistory = np.empty([(number_of_generations+1)*number_of_parents, number_of_parameters])
	#insert the value of initial parameters in history
	populationHistory[0:number_of_parents, :] = population

	for generation in range(number_of_generations):
		print("This is number %s generation" % (generation))
		if model_name=="RandomForestClassifier" and generation>50:
			break
		
		#train the dataset and obtain fitness
		#print("y_test",y_test)
		train_population_to_call = getattr(diff_models_train_genetic_opti, str(model_name)+"_train_pop")
		fitnessValue = train_population_to_call(population, X_train,X_test,y_train, y_test, list_encode_dics)
		print("fitnessValue")
		print(fitnessValue)
		#fitnessValue = train_population(population=population, dMatrixTrain=xgbDMatrix, dMatrixtest=xgbDMatrixTest, y_test=y_test)
		fitnessHistory[generation, :] = fitnessValue
		
		#best score in the current iteration
		print('Best MCC score in the this iteration = {}'.format(np.max(fitnessHistory[generation, :])))
		#survival of the fittest - take the top parents, based on the fitness value and number of parents needed to be selected
		parents = new_parents_selection(population=population, fitness=fitnessValue, numParents=number_of_parents_mating)
		
		#mate these parents to create children having parameters from these parents (we are using uniform crossover)
		children = crossover_uniform(parents=parents, childrenSize=(populationSize[0] - parents.shape[0], number_of_parameters))
		
		#add mutation to create genetic diversity
		mutate_to_call = getattr(diff_models_train_genetic_opti, str(model_name)+"_mutate")
		children_mutated = mutate_to_call(children, number_of_parameters)
		
		'''
		We will create new population, which will contain parents that where selected previously based on the
		fitness score and rest of them  will be children
		'''
		population[0:parents.shape[0], :] = parents #fittest parents
		population[parents.shape[0]:, :] = children_mutated #children
		
		populationHistory[(generation+1)*number_of_parents : (generation+1)*number_of_parents+ number_of_parents , :] = population #srore parent information


	#Best solution from the final iteration
	fitness = train_population_to_call(population, X_train,X_test,y_train, y_test, list_encode_dics)
	fitnessHistory[generation+1, :] = fitness#index of the best solution
	bestFitnessIndex = np.where(fitness == np.max(fitness))[0][0]#Best fitness

	print("Best fitness is =", fitness[bestFitnessIndex])#Best parameters
	param=population[bestFitnessIndex]
	params=zip(param_names, param)
	param_dictionary=dict(params)
	for key,value in param_dictionary.items():
		if key in [i[0] for i in list_encode_dics]:
			encode_param_index=[i[0] for i in list_encode_dics].index(key)
			(param_name,param_dic)=list_encode_dics[encode_param_index]
			param_dictionary[key]=param_dic[value]
	del fitnessHistory
	del populationHistory
	del population
	return fitness[bestFitnessIndex], param_dictionary, param_names

def get_param_dic(param_names, best_pop):
	param_dic={}
	for name,value in zip(param_names, best_pop):
		param_dic[name]=value

	return param_dic

def log_loss_func(weights,predictions, y_test):
   ''' scipy minimize will pass the weights as a numpy array '''
   final_prediction = 0
   for weight, prediction in zip(weights, predictions):
      #print(prediction)
      final_prediction += weight*prediction

   return log_loss(y_test, final_prediction)


#https://www.kaggle.com/hsperr/finding-ensamble-weights
#print("MODELS",models,"MODELS")
def optimize_weights(models,X_test, X_train, y_train, y_test):
	predictions = []
	for model_name, model in models:
		print("OPTIMIZE WEIGHT OF : ", model_name)
		try:
			model_params=model.get_params()
			model_params["probability"]=True
			model=model.set_params(**model_params)
		except:
			pass
		try:
			model_params=model.get_params()
			try:
				n_neighbors_param=model_params["n_neighbors"]
				model_params["n_neighbors"]=int(n_neighbors_param)
			except:
				pass
			try:
				n_estimators_param=model_params["n_estimators"]
				model_params["n_estimators"]=int(n_estimators_param)
			except:
				pass
			try:
				min_samples_split_param=model_params["min_samples_split"]
				model_params["min_samples_split"]=int(min_samples_split_param)
			except:
				pass
			model=model.set_params(**model_params)
		except:
			pass
		model.fit(X_train, np.ravel(y_train))
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

def get_voting_params(models,X_test, X_train, y_train, y_test):
   tup_dict = dict(models) 
   tup_dict.pop('Voting_all_soft')
   models = tuple(tup_dict.items())
   param_space=optimize_weights(models,X_test, X_train, y_train, y_test)
   weights= param_space
   model_dic={}
   model_dic['weights']=weights
   model_dic['estimators']=models

   return model_dic

def opti_run(models,X_train, y_train, X_test, y_test):
	model_dic={} #add all models here
	model_list=[]
	number_of_parameters_dic={"XGBoost":7,"LogisticRegression":7, "KNeighborsClassifier":4,"DecisionTreeClassifier":4, "GaussianNB":1,
	"RandomForestClassifier":7, "SVC":7}

	for model_name,model in models:
		model_list.append(model)

	for model_name, model in models:
		print("NOW OPTIMIZING: ",model_name)
		if model_name=="SVC" :  #model_name=="RandomForestClassifier" or
			index=models.index((model_name, model))
			param_dic={}
			continue
		if model_name=="Voting_all_soft":
			param_dic=get_voting_params(models,X_test, X_train, y_train, y_test)
			index=models.index((model_name, model))
		else:
			number_of_parameters=number_of_parameters_dic[model_name]
			index=models.index((model_name, model))
			try:
				best_fitness, param_dic, param_names=implement_opti(X_train,y_train,X_test,y_test, model, model_name, number_of_parameters)
				print(param_dic)
			except:
				index=models.index((model_name, model))
				param_dic={}

		
		#param_dic=get_param_dic(param_names, best_params)
		model=model.set_params(**param_dic)
		models[index]=(model_name,model)
	return models
		