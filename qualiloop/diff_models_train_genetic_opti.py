from sklearn.metrics import matthews_corrcoef
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import VotingClassifier
import xgboost as xgb
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef
from sklearn.exceptions import ConvergenceWarning, FitFailedWarning
from sklearn.utils._testing import ignore_warnings 
from sklearn.exceptions import NotFittedError
import gc



def fitness_mcc_score(y_true, y_pred):
	counter=0
	while counter<10:
		try:
			fitness=matthews_corrcoef(y_true, y_pred, sample_weight=None ) #maybe round((..),4)
			gc.collect()
			return fitness#train the data annd find fitness score
		except TypeError:
			print(TypeError)
			counter+=1
	return 0


def XGBoost_train_pop(population, X_train,X_test,y_train, y_test, list_encode_dics):
	X_train, X_test, y_train, y_test = X_train.to_numpy(), X_test.to_numpy(), y_train.to_numpy(), y_test.to_numpy()
	full_y=np.concatenate((np.unique(y_train), np.unique(y_test)))
	num_class=len(np.unique(full_y))
	
	dMatrixTrain = xgb.DMatrix(X_train, label=y_train)
	dMatrixTest = xgb.DMatrix(X_test, label=y_test)
	mcc_score = []
	print(np.unique(full_y))
	

	for i in range(population.shape[0]):
			gc.collect()
			param = { 'objective':'binary:logistic',
							'learning_rate': population[i][0],
							'n_estimators': population[i][1], 
							'max_depth': int(population[i][2]), 
							'min_child_weight': population[i][3],
							'gamma': population[i][4], 
							'subsample': population[i][5],
							'colsample_bytree': population[i][6],
							'seed': 24}
			if num_class>2:
				param['num_class']=num_class+1
				param['objective']='multi:softprob'
				print("num_class: ", num_class)
			num_round = 100
			xgb_classifier=xgb.XGBClassifier()

			try:
				with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
					xgbT = xgb.train(param, dMatrixTrain, num_round)
					preds = xgbT.predict(dMatrixTest)
					if num_class<=2:
						preds = preds>0.5  # transorms into true, false -> as xgboost predicts a probability by default
					else:
						nom_preds=[]
						for arr in preds:
							result = np.where(arr == np.amax(arr))
							nom_pred=result
							if len(nom_pred)>1:
								print("ERROR, nom_pred is greater than 1")
								nom_pred=nom_pred[0]
							nom_preds=np.append(nom_preds, nom_pred)
						preds=nom_preds
						tuple(int(x+1) for x in preds)

					print(tuple(preds))
					mcc_score.append(fitness_mcc_score(tuple(y_test), tuple(preds)))
			except NotFittedError:
				continue
			except ValueError:
				error_count+=1
			#	continue
			#except TypeError:
			#	continue
	return mcc_score
	


def XGBoost_mutate(crossover, number_of_parameters):
	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue[0:] = [0.01, 1.0] #min/max learning rate
	minMaxValue[1, :] = [10, 2000] #min/max n_estimator
	minMaxValue[2, :] = [1, 15] #min/max depth
	minMaxValue[3, :] = [0, 10.0] #min/max child_weight
	minMaxValue[4, :] = [0.01, 10.0] #min/max gamma
	minMaxValue[5, :] = [0.01, 1.0] #min/maxsubsample
	minMaxValue[6, :] = [0.01, 1.0] #min/maxcolsample_bytree
 
	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	if parameterSelect == 0: #learning_rate
		mutationValue = round(np.random.uniform(-0.5, 0.5), 2)
	if parameterSelect == 1: #n_estimators
		mutationValue = np.random.randint(-200, 200, 1)
	if parameterSelect == 2: #max_depth
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 3: #min_child_weight
		mutationValue = round(np.random.uniform(5, 5), 2)
	if parameterSelect == 4: #gamma
		mutationValue = round(np.random.uniform(-2, 2), 2)
	if parameterSelect == 5: #subsample
		mutationValue = round(np.random.uniform(-0.5, 0.5), 2)
	if parameterSelect == 6: #colsample
		mutationValue = round(np.random.uniform(-0.5, 0.5), 2)
	
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]    
	gc.collect()
	return crossover

def LogisticRegression_train_pop(population, X_train,X_test,y_train, y_test, list_encode_dics):
	mcc_score = []
	nr_for_loops=population.shape[0]
	error_count=0
	while nr_for_loops > 0 and error_count<10:
		for i in range(population.shape[0]):
			gc.collect()
			param = { #"penalty":population[i][0],
								"C":population[i][0],
								"fit_intercept":population[i][1], 
								"intercept_scaling":population[i][2], 
								"class_weight":population[i][3],
								"solver":population[i][4],
								"max_iter":int(population[i][5]), 
			"multi_class":population[i][6]}
			for key,value in param.items():
				if key in [i[0] for i in list_encode_dics]:
					encode_param_index=[i[0] for i in list_encode_dics].index(key)
					(param_name,param_dic)=list_encode_dics[encode_param_index]
					param[key]=param_dic[int(value)]
			model=LogisticRegression() 
			
			
			try:
				with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
					model=model.set_params(**param)
					model.fit(X_train, tuple(y_train))
					preds=model.predict(X_test)
					print(preds)
					score=fitness_mcc_score(y_test, tuple(preds))
					mcc_score.append(score)
					print("YEESSS")
					print("nr_for_loops", nr_for_loops)
					nr_for_loops=nr_for_loops-1
					if nr_for_loops <=0:
						break
			except NotFittedError:
				pass
			except ValueError:
				print(ValueError)
				error_count+=1
				pass
			#except TypeError:
			#	pass
			
	return mcc_score

def LogisticRegression_mutate(crossover, number_of_parameters):
	"""penalty : {'l1', 'l2', 'elasticnet', 'none'}, default='l2'
		dual : bool, default=False
		tol : float, default=1e-4
		C : float, default=1.0
		fit_intercept : bool, default=True
		intercept_scaling : float, default=1
		class_weight : dict or 'balanced', default=None
		random_state : int, RandomState instance, default=None
		solver : {'newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga'},             default='lbfgs'
		max_iter : int, default=100
		multi_class : {'auto', 'ovr', 'multinomial'}, default='auto'
		verbose : int, default=0
		warm_start : bool, default=False
		n_jobs : int, default=None
		l1_ratio : float, default=None"""

	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	#minMaxValue[0:] = [0, 3] #min/max penalty
	minMaxValue[0, :] = [0.1, 3] #min/max C
	minMaxValue[1, :] = [0, 1] #min/max fit_intercept
	minMaxValue[2, :] = [0.01, 10.0] #min/max intercept_scaling
	minMaxValue[3, :] = [0, 1] #min/max class_weight
	minMaxValue[4, :] = [0, 3] #min/max solver
	minMaxValue[5, :] = [20, 500] #min/max max_iter
	minMaxValue[6, :] = [0, 2] #min/max multi_class
 
	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	#if parameterSelect == 0:
	#	mutationValue = np.random.randint(-3, 3, 1)
	
	if parameterSelect == 0:
		mutationValue = round(np.random.uniform(-1,1), 2)
	if parameterSelect == 1: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 2: 
		mutationValue = round(np.random.uniform(-2,2), 2)
	if parameterSelect == 3: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 4: 
		mutationValue = np.random.randint(-2, 2, 1)
	if parameterSelect == 5: 
		mutationValue = np.random.randint(-100, 100, 1)
	if parameterSelect == 6: 
		mutationValue = np.random.randint(-2, 2, 1)
	
	
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]    
	gc.collect()
	return crossover

def KNeighborsClassifier_train_pop(population, X_train,X_test,y_train, y_test,list_encode_dics):
	mcc_score = []
	nr_for_loops=population.shape[0]
	error_count=0
	while nr_for_loops > 0 and error_count<10:
		for i in range(population.shape[0]):
				gc.collect()
				param = { "n_neighbors":int(population[i][0]), 
									"weights":population[i][1], 
									"algorithm":population[i][2], 
									"leaf_size":population[i][3]}
				for key,value in param.items():
					if key in [i[0] for i in list_encode_dics]:
						encode_param_index=[i[0] for i in list_encode_dics].index(key)
						(param_name,param_dic)=list_encode_dics[encode_param_index]
						param[key]=param_dic[value]
				model=KNeighborsClassifier()
				print(param)
				try:
					with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
						model=model.set_params(**param)
						model.fit(X_train, np.ravel(y_train))
						preds=model.predict(X_test)
						mcc_score.append(fitness_mcc_score(y_test, tuple(preds)))
						print("YESSS")
						nr_for_loops=nr_for_loops-1
						if nr_for_loops <=0:
							break
				except NotFittedError:
					continue
				except ValueError:
					print(ValueError)
					error_count+=1
				#except TypeError:
				#	continue
	return mcc_score

def KNeighborsClassifier_mutate(crossover, number_of_parameters):
	"""n_neighbors : int, default=5
		weights : {'uniform', 'distance'} or callable, default='uniform'
- 'uniform' : uniform weights.  All points in each neighborhood
- 'distance' : weight points by the inverse of their distance.
- [callable] : a user-defined function which accepts an
		algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, default='auto'
		leaf_size : int, default=30
		p : int, default=2
		metric : str or callable, default='minkowski'
		metric_params : dict, default=None
		n_jobs : int, default=None"""
	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue[0:] = [3, 20] #min/max n_neighbors
	minMaxValue[1, :] = [0, 1] #min/max weights
	minMaxValue[2, :] = [0, 3] #min/max algorithm
	minMaxValue[3, :] = [5, 100] #min/max leaf_size
	
	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	if parameterSelect == 0: 
		mutationValue = np.random.randint(-4, 4, 1)
	if parameterSelect == 1: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 2: 
		mutationValue = np.random.randint(-3, 3, 1)
	if parameterSelect == 3: 
		mutationValue = np.random.randint(-10, 10, 1)
 
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]    
	gc.collect()
	return crossover


def DecisionTreeClassifier_train_pop(population, X_train,X_test,y_train, y_test, list_encode_dics):
	mcc_score = []
	nr_for_loops=population.shape[0]
	error_count=0
	while nr_for_loops > 0 and error_count<10:
		for i in range(population.shape[0]):
				gc.collect()
				param = { "criterion":population[i][0], 
									"splitter":population[i][1], 
									"max_depth":int(population[i][2]), 
									"random_state":population[i][3]}
				for key,value in param.items():
					if key in [i[0] for i in list_encode_dics]:
						encode_param_index=[i[0] for i in list_encode_dics].index(key)
						(param_name,param_dic)=list_encode_dics[encode_param_index]
						param[key]=param_dic[value]
				model=DecisionTreeClassifier()
				print(param)
				try:
					with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
						model=model.set_params(**param)
						model.fit(X_train, np.ravel(y_train))
						preds=model.predict(X_test)
						print(tuple(preds))
						mcc_score.append(fitness_mcc_score(y_test, tuple(preds)))
						print("YESS")
						nr_for_loops=nr_for_loops-1
						if nr_for_loops <=0:
							break
				except NotFittedError:
					continue
				except ValueError:
					print(ValueError)
					error_count+=1
				except:
					continue
				#except TypeError:
				#	continue
	return mcc_score

def DecisionTreeClassifier_mutate(crossover, number_of_parameters):
	"""criterion : {"gini", "entropy"}, default="gini"
		splitter : {"best", "random"}, default="best"
		max_depth : int, default=None
		min_samples_split : int or float, default=2
		min_samples_leaf : int or float, default=1
		min_weight_fraction_leaf : float, default=0.0
		max_features : int, float or {"auto", "sqrt", "log2"}, default=None
		random_state : int, RandomState instance or None, default=None
		max_leaf_nodes : int, default=None
		min_impurity_decrease : float, default=0.0
		min_impurity_split : float, default=0
		class_weight : dict, list of dict or "balanced", default=None
	ccp_alpha : non-negative float, default=0.0"""
	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue[0:] = [0, 1] #min/max criterion
	minMaxValue[1, :] = [0, 1] #min/max splitter
	minMaxValue[2, :] = [1, 30] #min/max max_depth
	minMaxValue[3, :] = [0, 1] #min/max random_state
 
	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	if parameterSelect == 0:
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 1: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 2: 
		mutationValue = np.random.randint(-5, 5, 1)
	if parameterSelect == 3: 
		mutationValue = np.random.randint(-1, 1, 1)
	
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]    
	gc.collect()
	return crossover


def GaussianNB_train_pop(population, X_train,X_test,y_train, y_test, list_encode_dics):
	mcc_score = []
	nr_for_loops=population.shape[0]
	error_count=0
	while nr_for_loops > 0 and error_count<10:
		for i in range(population.shape[0]):
				gc.collect()
				param = { "var_smoothing":population[i][0]}
				model=GaussianNB()
				print(param)
				try:
					with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
						model=model.set_params(**param)
						model.fit(X_train, np.ravel(y_train))
						preds=model.predict(X_test)
						print(y_test, tuple(preds))
						mcc_score.append(fitness_mcc_score(y_test, tuple(preds)))
						print("YESS")
						nr_for_loops=nr_for_loops-1
						if nr_for_loops <=0:
							break
				except NotFittedError:
					continue
				except ValueError:
					print(ValueError)
					error_count+=1
				#except TypeError:
				#	continue
	return mcc_score

def GaussianNB_mutate(crossover, number_of_parameters):
	"""GaussianNB
		priors : array-like of shape (n_classes,)
		var_smoothing : float, default=1e-9"""

	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue[0:] = [float(1e-11), float(1e-2)] #min/max var_smoothing
	
	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	if parameterSelect == 0:
		mutationValue = np.random.uniform(float(-1e-4), float(1e-4))
	
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]    
	gc.collect()
	return crossover

def RandomForestClassifier_train_pop(population, X_train,X_test,y_train, y_test, list_encode_dics):
	mcc_score = []
	nr_for_loops=population.shape[0]
	error_count=0
	while nr_for_loops > 0 and error_count<10:
		for i in range(population.shape[0]):
				gc.collect()
				param = { "n_estimators":int(population[i][0]), 
									"criterion":population[i][1], 
									"min_samples_split":int(population[i][2]), 
									"bootstrap":population[i][3], 
									"random_state":population[i][4], 
									"class_weight":population[i][5],
									"max_depth": int(population[i][6]), 
				}
				for key,value in param.items():
					if key in [i[0] for i in list_encode_dics]:
						encode_param_index=[i[0] for i in list_encode_dics].index(key)
						(param_name,param_dic)=list_encode_dics[encode_param_index]
						param[key]=param_dic[value]
				model=RandomForestClassifier()
				print(param)
				try:
					with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
						model=model.set_params(**param)
						model.fit(X_train, tuple(y_train))
						preds=model.predict(X_test)
						print(y_test, tuple(preds))
						try:
							mcc_score.append(fitness_mcc_score(tuple(y_test.to_numpy()), tuple(preds)))
							print("YESSS")
							nr_for_loops=nr_for_loops-1
						except:
							break
						if nr_for_loops <=0:
							break
				except NotFittedError:
					continue
				except ValueError:
					print(ValueError)
					error_count+=1
				#except TypeError:
				#	continue
	return mcc_score
	
def RandomForestClassifier_mutate(crossover, number_of_parameters):
	"""n_estimators : int, default=100
		criterion : {"gini", "entropy"}, default="gini"
		max_depth : int, default=None
		min_samples_split : int or float, default=2
		min_samples_leaf : int or float, default=1
		min_weight_fraction_leaf : float, default=0.0
		max_features : {"auto", "sqrt", "log2"}, int or float, default="auto"
		max_leaf_nodes : int, default=None
		min_impurity_decrease : float, default=0.0
		min_impurity_split : float, default=None
		bootstrap : bool, default=True
		oob_score : bool, default=False
		n_jobs : int, default=None
		random_state : int, RandomState instance or None, default=None
		verbose : int, default=0
		warm_start : bool, default=False
		class_weight : {"balanced", "balanced_subsample"}, dict or list of dicts,             default=None
		ccp_alpha : non-negative float, default=0.0
		max_samples : int or float, default=None"""

	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue[0:] = [50, 200] #min/max n_estimators
	minMaxValue[1, :] = [0, 1] #min/max criterion
	minMaxValue[2, :] = [2, 10] #min/max min_samples_split
	minMaxValue[3, :] = [0, 1] #min/max bootstrap
	minMaxValue[4, :] = [0, 1] #min/max random_state
	minMaxValue[5, :] = [0, 2] #min/max class_weight
	minMaxValue[6, :] = [1, 20] #min/max max_depth

	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	if parameterSelect == 0:
		mutationValue = np.random.randint(-50, 50, 1)
	if parameterSelect == 1: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 2:
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 3: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 4: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 5: 
		mutationValue = np.random.randint(-2, 2, 1)
	if parameterSelect == 6: 
		mutationValue = np.random.randint(-5, 5, 1)
	
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]    
	gc.collect()
	return crossover

def SVC_train_pop(population, X_train,X_test,y_train, y_test, list_encode_dics):
	mcc_score = []
	nr_for_loops=population.shape[0]
	error_count=10
	while nr_for_loops > 0 and error_count<10:
		for i in range(population.shape[0]):
				gc.collect()
				param = { "C":population[i][0], 
									"kernel":population[i][1], 
									"degree":population[i][2], 
									"gamma":population[i][3], 
									"coef0":population[i][4] ,  
									"class_weight": population[i][5], 
									"random_state":population[i][6],
									"max_iter":int(-1),
									"probability":True,
									"decision_function_shape": "ovo"}
				for key,value in param.items():
					if key in [i[0] for i in list_encode_dics]:
						encode_param_index=[i[0] for i in list_encode_dics].index(key)
						(param_name,param_dic)=list_encode_dics[encode_param_index]
						param[key]=param_dic[value]
				model=SVC()
				print(param)
				try:
					with ignore_warnings(category=[ConvergenceWarning, FitFailedWarning]):
						model=model.set_params(**param)
						model.fit(X_train, np.ravel(y_train))
						preds=model.predict(X_test)
						try:
							mcc_score.append(fitness_mcc_score(tuple([int(i) for i in y_test.to_numpy()]), tuple([int(i) for i in preds])))
							print("YESS")
							nr_for_loops=nr_for_loops-1
						except:
							continue
						if nr_for_loops <=0:
							break
				except NotFittedError:
					continue
				except ValueError:
					print(ValueError)
					error_count+=1
				#except TypeError:
				#	continue
	return mcc_score

def SVC_mutate(crossover, number_of_parameters):
	""" C : float, default=1.0
		kernel : {'linear', 'poly', 'rbf', 'sigmoid', 'precomputed'}, default='rbf'
		degree : int, default=3
		gamma : {'scale', 'auto'} or float, default='scale'
		coef0 : float, default=0.0
		shrinking : bool, default=True
		probability : bool, default=False
		tol : float, default=1e-3
		cache_size : float, default=200
		class_weight : dict or 'balanced', default=None
		verbose : bool, default=False
		max_iter : int, default=-1
		decision_function_shape : {'ovo', 'ovr'}, default='ovr'
		break_ties : bool, default=False
	random_state : int, RandomState instance or None, default=None"""
	#Define minimum and maximum values allowed for each parameterminMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue = np.zeros((number_of_parameters, 2))
	minMaxValue[0:] = [0.01, 1.0] #min/max C
	minMaxValue[1, :] = [0, 3] #min/max kernel
	minMaxValue[2, :] = [0, 10] #min/max degree
	minMaxValue[3, :] = [0, 1] #min/max gamma
	minMaxValue[4, :] = [0.01, 10.0] #min/max coef0
	minMaxValue[5, :] = [0, 1] #min/max class_weight
	minMaxValue[6, :] = [0, 1] #min/max random_state
	#changes a single gene in each offspring randomly.
	mutationValue = 0
	parameterSelect = np.random.randint(0, number_of_parameters, 1)
	print(parameterSelect)
	if parameterSelect == 0:
		mutationValue = round(np.random.uniform(-0.2, 0.2), 2)
	if parameterSelect == 1: 
		mutationValue = np.random.randint(-3, 3, 1)
	if parameterSelect == 2:
		mutationValue = np.random.randint(-5, 5, 1)
	if parameterSelect == 3: 
		mutationValue = round(np.random.randint(-1, 1), 2)
	if parameterSelect == 4: 
		mutationValue = round(np.random.uniform(-5, 5), 2)
	if parameterSelect == 5: 
		mutationValue = np.random.randint(-1, 1, 1)
	if parameterSelect == 6: 
		mutationValue = np.random.randint(-1, 1, 1)
	
	#indtroduce mutation by changing one parameter, and set to max or min if it goes out of range
	for idx in range(crossover.shape[0]):
		crossover[idx, parameterSelect] = crossover[idx, parameterSelect] + mutationValue
		if(crossover[idx, parameterSelect] > minMaxValue[parameterSelect, 1]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 1]
		if(crossover[idx, parameterSelect] < minMaxValue[parameterSelect, 0]):
			crossover[idx, parameterSelect] = minMaxValue[parameterSelect, 0]
	gc.collect()    
	return crossover