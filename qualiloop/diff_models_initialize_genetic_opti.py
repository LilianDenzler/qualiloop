import sys
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./CDRH3lib')
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef
import xgboost as xgb
import random


def encode_cat(param_name, val_list):
    #encode categorical parameters as integers,
    #gives a dictionary for decoding
    param_dic={}
    param_len=len(val_list)
    for i, val in zip(range(0,param_len),val_list):
        param_dic[i]=val
    return param_len, param_dic

def XGBoost_initialize_pop(numberOfParents):
    learningRate = np.empty([numberOfParents, 1])
    nEstimators = np.empty([numberOfParents, 1], dtype = np.uint8)
    maxDepth = np.empty([numberOfParents, 1], dtype = np.uint8)
    minChildWeight = np.empty([numberOfParents, 1])
    gammaValue = np.empty([numberOfParents, 1])
    subSample = np.empty([numberOfParents, 1])
    colSampleByTree =  np.empty([numberOfParents, 1])
    for i in range(numberOfParents):
        #print(i)
        learningRate[i] = round(random.uniform(0.01, 1), 2)
        nEstimators[i] = random.randrange(10, 1500, step = 25)
        maxDepth[i] = int(random.randrange(1, 10, step= 1))
        minChildWeight[i] = round(random.uniform(0.01, 10.0), 2)
        gammaValue[i] = round(random.uniform(0.01, 10.0), 2)
        subSample[i] = round(random.uniform(0.01, 1.0), 2)
        colSampleByTree[i] = round(random.uniform(0.01, 1.0), 2)
    
    population = np.concatenate((learningRate, nEstimators, maxDepth, minChildWeight, gammaValue, subSample, colSampleByTree), axis= 1)
    param_names=["learningRate", "nEstimators", "maxDepth", "minChildWeight", "gammaValue", "subSample","colSampleByTree"]
    list_encode_dics=[]
    return population, param_names, list_encode_dics

def LogisticRegression_initialize_pop(numberOfParents):
    #penalty=np.empty([numberOfParents, 1], dtype = np.uint8)
    #param_len_penalty, param_dic_penalty=encode_cat("penalty", ['l1', 'l2', 'elasticnet', 'none'])

    C=np.empty([numberOfParents, 1])

    fit_intercept=np.empty([numberOfParents, 1], dtype = np.uint8)
    param_len_fit_intercept, param_dic_fit_intercept=encode_cat("fit_intercept", [True, False])

    intercept_scaling=np.empty([numberOfParents, 1])

    class_weight=np.empty([numberOfParents, 1], dtype = np.uint8)
    param_len_class_weight, param_dic_class_weight=encode_cat("class_weight", ['balanced',None])

    solver=np.empty([numberOfParents, 1], dtype = np.uint8)
    param_len_solver, param_dic_solver=encode_cat("solver", ['newton-cg', 'lbfgs', 'sag', 'saga'])

    max_iter=np.empty([numberOfParents, 1])

    multi_class=np.empty([numberOfParents, 1], dtype = np.uint8)
    param_len_multi_class, param_dic_multi_class=encode_cat("multi_class", ['auto', 'ovr', 'multinomial'])

    for i in range(numberOfParents):
        #penalty[i] = random.randrange(0, (param_len_penalty))
        C[i] = round(random.uniform(0, 5), 2)
        fit_intercept[i]=random.randrange(0, (param_len_fit_intercept))
        intercept_scaling[i] = round(random.uniform(0, 10), 2)
        class_weight[i]=random.randrange(0, (param_len_class_weight))
        solver[i]=random.randrange(0, (param_len_solver))
        max_iter[i] = int(random.randrange(1, 500, step= 1))
        multi_class[i]=random.randrange(0, (param_len_multi_class))
    population = np.concatenate(( C, fit_intercept, intercept_scaling, class_weight, solver,max_iter, multi_class), axis= 1) #penalty
    param_names=["C", "fit_intercept", "intercept_scaling", "class_weight","solver","max_iter", "multi_class"]
    list_encode_dics=[ ("fit_intercept", param_dic_fit_intercept), 
                      ("class_weight", param_dic_class_weight), ("solver", param_dic_solver), ("multi_class", param_dic_multi_class)]#("penalty",param_dic_penalty)
    return population, param_names, list_encode_dics

def  LinearDiscriminantAnalysis_initialize_pop(numberOfParents):
    solver=np.empty([numberOfParents, 1], dtype = np.uint8)
    param_len_solver, param_dic_solver=encode_cat("solver", ['svd', 'lsqr', 'eigen'])
    for i in range(numberOfParents):
        solver[i] = random.randrange(1, (param_len_solver+1))
    population=np.concatenate((solver), axis=1)
    param_names=["solver"]
    list_encode_dics=[("solver", param_dic_solver)]
    return population,param_names, list_encode_dics
    
    #'shrinkage':['auto',dists.uniform(0, 1)],
    #'tol': tuple(np.linspace(0, 500, 100))


def KNeighborsClassifier_initialize_pop(numberOfParents):
    n_neighbors=np.empty([numberOfParents, 1])

    weights=np.empty([numberOfParents, 1])
    param_len_weights, param_dic_weights=encode_cat("weights", ['uniform', 'distance'])

    algorithm=np.empty([numberOfParents, 1])
    param_len_algorithm, param_dic_algorithm=encode_cat("algorithm", ['auto', 'ball_tree', 'kd_tree', 'brute'])

    leaf_size=np.empty([numberOfParents, 1])

    for i in range(numberOfParents):
        n_neighbors[i] = int(random.randrange(1, 20, step= 1))
        weights[i]=random.randrange(0, (param_len_weights))
        algorithm[i]=random.randrange(0, (param_len_algorithm))
        leaf_size[i] = int(random.randrange(1, 100, step= 1))

    population = np.concatenate((n_neighbors, weights, algorithm, leaf_size), axis= 1)
    param_names=["n_neighbors", "weights", "algorithm", "leaf_size"]
    list_encode_dics=[("weights", param_dic_weights), ("algorithm", param_dic_algorithm)]
    return population, param_names, list_encode_dics

def DecisionTreeClassifier_initialize_pop(numberOfParents):
    criterion=np.empty([numberOfParents, 1])
    param_len_criterion, param_dic_criterion=encode_cat("criterion", ['gini', 'entropy'])

    splitter=np.empty([numberOfParents, 1])
    param_len_splitter, param_dic_splitter=encode_cat("splitter", ['best', 'random'])

    max_depth=np.empty([numberOfParents, 1])

    random_state=np.empty([numberOfParents, 1])
    param_len_random_state, param_dic_random_state=encode_cat("random_state", [None, 42])

    for i in range(numberOfParents):
        criterion[i]=random.randrange(0, (param_len_criterion))
        splitter[i]=random.randrange(0, (param_len_criterion))
        max_depth[i] = int(random.randrange(1, 30))
        random_state[i]=random.randrange(0, (param_len_criterion))
    #'max_depth': [dists.randint(0,500),None],
    #'min_samples_leaf': dists.uniform(0, 0.5),
    #'max_features': ['auto', 'sqrt', 'log2',dists.uniform(0, 1)],
    #'max_leaf_nodes': [dists.randint(0,500),None],
    #'min_impurity_decrease': dists.uniform(0, 500),
    #'min_impurity_split': dists.uniform(0, 500),
    #'ccp_alpha': dists.uniform(0, 500),
    population = np.concatenate((criterion,splitter, max_depth, random_state), axis= 1)
    param_names=["criterion", "splitter", "max_depth", "random_state"]
    list_encode_dics=[("criterion", param_dic_criterion), ("splitter", param_dic_splitter), ("random_state", param_dic_random_state)] 
    return population, param_names, list_encode_dics

def GaussianNB_initialize_pop(numberOfParents):
    var_smoothing=np.empty([numberOfParents, 1])
    for i in range(numberOfParents):
        var_smoothing[i] = round(random.uniform(0, 0.5), 2)
    population=var_smoothing
    param_names=["var_smoothing"]
    list_encode_dics=list()
    return population, param_names, list_encode_dics
    

def RandomForestClassifier_initialize_pop(numberOfParents):
    n_estimators=np.empty([numberOfParents, 1])

    criterion=np.empty([numberOfParents, 1])
    param_len_criterion, param_dic_criterion=encode_cat("criterion", ['gini', 'entropy'])

    min_samples_split=np.empty([numberOfParents, 1])

    bootstrap=np.empty([numberOfParents, 1])
    param_len_bootstrap, param_dic_bootstrap=encode_cat("bootstrap", [True,False])

    random_state=np.empty([numberOfParents, 1])
    param_len_random_state, param_dic_random_state=encode_cat("random_state", [None,42])

    class_weight=np.empty([numberOfParents, 1])
    param_len_class_weight, param_dic_class_weight=encode_cat("class_weight", ['balanced', 'balanced_subsample', None])

    max_depth=np.empty([numberOfParents, 1])
    for i in range(numberOfParents):
        n_estimators[i] = int(random.randrange(50, 200, step= 1))
        criterion[i] = random.randrange(0, (param_len_criterion))
        min_samples_split[i] = int(random.randrange(2, 10))
        bootstrap[i] = random.randrange(0, (param_len_bootstrap))
        random_state[i] = random.randrange(0, (param_len_random_state))
        class_weight[i] = random.randrange(0, (param_len_class_weight))
        max_depth[i] = int(random.randrange(1, 20))

    population=np.concatenate((n_estimators, criterion, min_samples_split, bootstrap, random_state, class_weight, max_depth), axis=1)
    param_names=["n_estimators", "criterion", "min_samples_split", "bootstrap", "random_state", "class_weight", "max_depth"]
    list_encode_dics=[("criterion", param_dic_criterion), ("bootstrap", param_dic_bootstrap), ("random_state", param_dic_random_state), ("class_weight", param_dic_class_weight)]
    return population, param_names, list_encode_dics
      #'max_depth': [dists.randint(0,500),None],
      #'min_samples_leaf': dists.uniform(0, 0.5),
      #'min_weight_fraction_leaf': dists.randint(0,0.5),
      #'max_features':['auto', 'sqrt', 'log2',dists.uniform(0, 1)],
      #'max_leaf_nodes': [dists.randint(0,500),None],
      #'min_impurity_decrease': dists.uniform(0, 500),
      #'min_impurity_split': dists.uniform(0, 500),
      #'ccp_alpha': dists.uniform(0, 500),
      #'max_samples': [dists.uniform(0, 1),None]

def SVC_initialize_pop(numberOfParents):
    C=np.empty([numberOfParents, 1])

    kernel=np.empty([numberOfParents, 1])
    param_len_kernel, param_dic_kernel=encode_cat("kernel", ['linear', 'poly', 'rbf', 'sigmoid'])

    degree=np.empty([numberOfParents, 1])

    gamma=np.empty([numberOfParents, 1])
    param_len_gamma, param_dic_gamma=encode_cat("gamma", ['scale', 'auto'])

    coef0=np.empty([numberOfParents, 1])

    class_weight=np.empty([numberOfParents, 1])
    param_len_class_weight, param_dic_class_weight=encode_cat("class_weight", ['balanced', None])

    random_state=np.empty([numberOfParents, 1])
    param_len_random_state, param_dic_random_state=encode_cat("random_state", [None,42])
    for i in range(numberOfParents):
        C[i] = round(random.uniform(0, 1), 2)
        kernel[i] = random.randrange(0, (param_len_kernel))
        degree[i] = int(random.randrange(0, 15))
        gamma[i] = random.randrange(0, (param_len_gamma))
        coef0[i]= round(random.uniform(0, 10), 2)
        class_weight[i] = random.randrange(0, (param_len_class_weight))
        random_state[i] = random.randrange(0, (param_len_random_state))
    population=np.concatenate((C,kernel,degree, gamma, coef0, class_weight, random_state), axis=1)
    param_names=["C", "kernel", "degree", "gamma", "coef0", "class_weight", "random_state"]
    list_encode_dics=[("kernel", param_dic_kernel), ("gamma", param_dic_gamma), ("class_weight",param_dic_class_weight), ("random_state", param_dic_random_state)]
    return population, param_names, list_encode_dics
   
      #'shrinking': [True, False],
      #'tol': tuple(np.linspace(0, 500, 100)),
      #'break_ties': [True,False],