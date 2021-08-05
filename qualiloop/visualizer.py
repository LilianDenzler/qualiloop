#!/usr/bin/python3
"""/*************************************************************************
	 Program:    Qualiloop Visualizer
	 File:       visualizer.py
	 
	 Version:    V1.1
	 Date:       19.03.21
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
	 The program contains all functions needed for different forms of data
	 preprocessing and feature encoding.
	 All of these differnet pipelines will be used to assess which types of 
	 preprocessing works best with which model. 


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

from yellowbrick.features import Rank2D
from yellowbrick.features import RadViz
from yellowbrick.features.pcoords import parallel_coordinates
from yellowbrick.features import PCA
from yellowbrick.features import JointPlotVisualizer
from yellowbrick.target import ClassBalance
from yellowbrick.model_selection import ValidationCurve
from yellowbrick.model_selection import LearningCurve
from yellowbrick.target import FeatureCorrelation
from yellowbrick.model_selection import FeatureImportances

from yellowbrick.model_selection import RFECV

from yellowbrick.classifier import ROCAUC
from yellowbrick.classifier import PrecisionRecallCurve
from yellowbrick.classifier import ClassificationReport
from yellowbrick.classifier import ClassPredictionError
from yellowbrick.classifier import DiscriminationThreshold



import math
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import sys
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import matthews_corrcoef
from sklearn.preprocessing import MinMaxScaler
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import make_scorer

from sklearn.tree import DecisionTreeClassifier

from sklearn.ensemble import RandomForestClassifier


#plots a heatmap of where data is missing. DO before cleaning to show where nans are
#https://kaiserm.medium.com/how-to-tackle-multicollinearity-79afe58e9479
# simple heat map showing where we are missing data

###################################################################################################
#TOPLAYER
#########


def run(full_dataset_df, RMSD_mode):
	del_RMSD_cols=[col for col in RMSD_mode if "nom" not in col]
	target_y=[col for col in RMSD_mode if "nom" in col][0]
	full_dataset_df_nom=full_dataset_df[full_dataset_df.columns[~full_dataset_df.columns.isin(del_RMSD_cols)]]
	print(full_dataset_df.columns)
	print(full_dataset_df_nom.columns)
	correlation_plot_all(full_dataset_df_nom)
	y=full_dataset_df[str(target_y)]
	X=full_dataset_df_nom.drop(columns=target_y, errors='ignore')
	pearson_ranker(X,y)
	random_forest_features(X, y)


def get_X_y(full_dataset_df, RMSD_mode):
	features_list=[]
	for i in full_dataset_df.columns:
		if "nom" in i and i in RMSD_mode:
			target_y=i
		if "bin" in i or "local_AA" in i or "global_AA" in i or "local_CA" in i or "global_CA" in i:
			pass
		else:
			features_list.append(i)
	X=full_dataset_df[features_list].astype(float)
	y=full_dataset_df[target_y].astype(float)
	return X, y

###################################################################################################
#PREANALYSIS
#########

def class_distribution(full_dataset_df, RMSD_mode,bin_thresholds, save_name, run_dir=None):
	#RMSD_mode is list of all target y's of the RMSD mode, i.e. ["local_CA_bin1", ..."local_CA_nom"]
	dim_a=math.floor(math.sqrt(len(RMSD_mode)))
	dim_b=math.floor((len(RMSD_mode))/dim_a)
	fig, axes = plt.subplots(dim_a, dim_b)
	fig.tight_layout(h_pad=2, w_pad=3, pad=3)
	
	grid_plot=[]
	counter=0
	for b in range(dim_b):
		for a in range(dim_a):
			target_y=RMSD_mode[counter]
			print("DIMESNIONS",a,b)
			print("TARGET")
			
			if "bin" in target_y:
				print("target_y",target_y)
				class_labels=full_dataset_df[target_y].unique().tolist()
				plot = ClassBalance(labels=class_labels, ax=axes[a][b], size=(1620, 1080))
				y=full_dataset_df[target_y]
				plot.fit(y) 
				plot.finalize()
				axes[a][b].set_ylabel('Number of Instances')
				#axes[a][b].set_xlabel('Class Label')
				axes[a][b].set_title(target_y+' Class Distribution')
				counter+=1
				
			elif "nom" in target_y:
				y=full_dataset_df[target_y]
				class_labels=full_dataset_df[target_y].unique().tolist()
				plot = ClassBalance(labels=class_labels, ax=axes[a][b], size=(1620, 1080))
				plot.fit(y)
				plot.finalize()
				axes[a][b].set_ylabel('Number of Instances')
				#axes[a][b].set_xlabel('Class Label')
				axes[a][b].set_title(target_y+' Class Distribution')
				#axes[a][b].vlines(bin_thresholds, 0, 1000, linestyles='solid')
				counter+=1
				
	fig.suptitle('Class Distributions')
	plt.subplots_adjust(top=0.85)
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()

def num_distribution(full_dataset_df, RMSD_mode,bin_thresholds, val_dic,RMSD_name,save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	target_y=RMSD_name
	sns.set(font_scale=1)
	fig=sns.distplot(full_dataset_df[target_y])
	plt.xlabel(str(target_y)+ " RMSD in Angstrom")
	plt.xlim(0,None)
	plt.ylabel("Nr. of Models")
	plt.title("AbDb Generated Model Quality Distribution")
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	#plt.clf()

	nom_values=val_dic["nom_label"]
	threshlows=val_dic["lower_threshold"]
	threshups=val_dic["upper_threshold"]
	print("THRESHLOW")
	print(threshlows)
	print(threshups)
	colours=sns.color_palette("colorblind")
	counter=0
	for i in bin_thresholds:
		plt.axvline(i, color="k", linestyle="--", label="_"*counter+"Bin Thresholds")  #use counter as labels with _ are ignored -> so only one label for bins on graphs
		counter+=1
	print(nom_values.tolist(), threshlows.tolist(), threshups.tolist(),colours)
	for label,down,up, colour in zip(nom_values.tolist(), threshlows.tolist(), threshups.tolist(),colours):
		if up==100:
			up=max(full_dataset_df[target_y].to_numpy())
		print(down, up, label)
		plt.axvspan(down, up, alpha=0.3, color=colour, label="Nominal Class"+str(label))
	fig.legend()
	save_name=save_name.replace(".png", "")+"_nom_bin.png"
	plt.xticks(rotation=45)

	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()


#CLEANING PROCESS
def data_preanalysis_nans(full_dataset_df, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	full_dataset_df=full_dataset_df[[col for col in full_dataset_df.columns if "local_AA" not in col and "local_CA" not in col and "global_AA" not in col and "global_CA" not in col]]
	mask=~full_dataset_df.isin([np.nan, np.inf, -np.inf, "nan", "None"])
	heat_map = sns.heatmap(mask, yticklabels = False, cbar = True, cmap = "PuRd", vmin = 0, vmax = 1)
	plt.savefig(os.path.join(run_dir,"graphs",save_name+"preanalysis_nans.png"))
	plt.clf()

def feature_box_plot(full_dataset_df, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	full_dataset_df=full_dataset_df[[col for col in full_dataset_df.columns if "local_AA" not in col and "local_CA" not in col and "global_AA" not in col and "global_CA" not in col]]
	features=list(full_dataset_df.select_dtypes(include=[np.number]).columns.values)
	X=full_dataset_df[features]
	labels=X.columns
	fig = plt.figure(figsize=(15,7))
	plt.boxplot(X.to_numpy(), vert=True, patch_artist=True, labels=labels) 
	plt.ylabel('observed value')
	plt.xticks(rotation=45)
	plt.title('Box Plot of all Features')
	print(run_dir)
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()

###################################################################################################
#FEATURE ANALYSIS
#################

def single_feature_corr(full_dataset_df, RMSD_mode, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	dim_b=len(RMSD_mode)
	features_list=[]
	for i in full_dataset_df.columns:
		if "nom" in i or "bin" in i or "local_AA" in i or "global_AA" in i or "local_CA" in i or "global_CA" in i:
			pass
		else:
			features_list.extend([i])
	targets_list=[col for col in full_dataset_df.columns if col not in features_list]
	if not os.path.exists(os.path.join(run_dir,"graphs", "single_feature_corr/")):
		os.makedirs(os.path.join(run_dir,"graphs", "single_feature_corr/"))
		visualizer_out_path_sub=os.path.join(run_dir,"graphs", "single_feature_corr/")
	else:
		visualizer_out_path_sub=os.path.join(run_dir,"graphs", "single_feature_corr/")
	for target_y in targets_list:
		for feature in features_list:
			#if "tip" in feature or "blosum" in feature or "one_hot" in feature or "NLF" in feature:
			#	continue
			print(feature)
			plot=JointPlotVisualizer(columns=feature, size=(1620, 1080))
			X=full_dataset_df[features_list].astype(float)
			y=full_dataset_df[target_y].astype(float)
			plot.fit_transform(X, y)
			plot.finalize()
			plt.ylabel(target_y)
			plt.xlabel(feature)
			plt.title("Correlation of "+feature+" and"+target_y)
	
			plt.xticks(rotation=45)
			plt.savefig(os.path.join(visualizer_out_path_sub,target_y+feature+save_name))
			plt.clf()

"""def feature_corr_alt():
	#alternative to single_feature_corr
	X, y = data['data'], data['target']
	X_pd = pd.DataFrame(X, columns=data['feature_names'])

	# Create a list of the features to plot
	features = ['alcohol', 'ash', 'hue', 'proline', 'total_phenols']

	# Instaniate the visualizer
	visualizer = FeatureCorrelation(
			method='mutual_info-classification', feature_names=features, sort=True
	)"""


def rad_viz_analysis(full_dataset_df, RMSD_mode, save_name, run_dir=None):
	#https://www.scikit-yb.org/en/latest/api/features/radviz.html
	#visualize each feature as dimension on circle, diff colours for each class-> see if diff coloured points are separatable in circle
	#gives measure of how noisy dataset is, and how well instances are classified with given features
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	X,y=get_X_y(full_dataset_df, RMSD_mode)
	for i in RMSD_mode:
		if "nom" in i:
			target_y=i
	class_labels=full_dataset_df[target_y].unique().tolist()
	y=pd.Series(y)
	X=pd.DataFrame(X)
	X.reset_index(drop=True, inplace=True)
	visualizer = RadViz(classes=class_labels, alpha=0.5, size=(1620, 1080))
	visualizer.fit_transform(X, y) 
	visualizer.finalize()               
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()    


def pearson_ranker(full_dataset_df, RMSD_mode, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	X,y=get_X_y(full_dataset_df, RMSD_mode)
	visualizer = Rank2D(algorithm='pearson', size=(1620, 1080))
	visualizer.fit_transform(X, y)  
	visualizer.finalize()      
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()             


def covariance_ranker(full_dataset_df, RMSD_mode, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	X,y=get_X_y(full_dataset_df, RMSD_mode)
	visualizer = Rank2D(algorithm='covariance', size=(1620, 1080))
	visualizer.fit(X, y)          
	visualizer.transform(X) 
	visualizer.finalize()      
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()                


def correlation_plot_all(full_dataset_df, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	plt.figure()
	cor = full_dataset_df.corr()
	sns.heatmap(cor, cmap=plt.cm.Reds) #annot=True
	plt.savefig(os.path.join(run_dir,"graphs", save_name))
	plt.clf()  


def parallel_coords(full_dataset_df, RMSD_mode, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	X,y=get_X_y(full_dataset_df, RMSD_mode)
	for i in RMSD_mode:
		if "nom" in i:
			target_y=i
	class_labels=full_dataset_df[target_y].unique().tolist()
	visualizer = parallel_coordinates(X, y, classes=class_labels, features=X.columns)
	visualizer.finalize()
	plt.savefig(os.path.join(run_dir,"graphs",save_name)) 
	plt.clf()  

def PCA_decomposition(full_dataset_df, RMSD_mode, save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	from sklearn import datasets
	iris = datasets.load_iris()
	X = iris.data
	y = iris.target
	print(X)
	print(y)
	print(type(X), type(y))
	for i in RMSD_mode:
		if "nom" in i:
			target_y=i
	full_dataset_df=full_dataset_df.sort_values(target_y)
	features_list=[]
	for i in full_dataset_df.columns:
		if "nom" in i and i in RMSD_mode:
			target_y=i
		if "bin" in i or "local_AA" in i or "global_AA" in i or "local_CA" in i or "global_CA" in i:
			pass
		else:
			features_list.append(i)
	X=full_dataset_df[features_list].astype(float).to_numpy()
	y=full_dataset_df[target_y].astype(int).to_numpy()
	class_labels=full_dataset_df[target_y].unique().tolist()
	#class_labels=np.array([int(i) for i in class_labels])
	print(class_labels)
	print(RMSD_mode)
	print(y)
	visualizer = PCA(scale=True, proj_features=True, classes=class_labels, size=(1620, 1080))
	visualizer.fit(X,y)
	visualizer.transform(X)
	visualizer.finalize()
	plt.savefig(os.path.join(run_dir,"graphs",save_name)) 
	plt.clf()  
	visualizer = PCA(scale=True, proj_features=True, projection=3, classes=class_labels, size=(1620, 1080))
	visualizer.fit_transform(X, y)
	visualizer.finalize()
	plt.savefig(os.path.join(run_dir,"graphs","3D_"+save_name)) 
	plt.clf()  

####################################################################################################################
#MODEL TRAINING
###############
def ROCAUC_curve(X_train, y_train, X_test, y_test, RMSD_mode, model,save_name,run_name=None, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs", run_name)):
		os.makedirs(os.path.join(run_dir,"graphs", run_name))
	visualizer_out_path_sub=os.path.join(run_dir,"graphs", run_name)
	if isinstance(y_train, (np.ndarray)) ==False:
		y_train=y_train.to_numpy()
	class_labels=np.unique(y_train).tolist()
	visualizer = ROCAUC(model, classes=class_labels)
	print(model)
	try:
		param={}
		param['probability']=False
		model=model.set_params(**param)
	except:
		pass
	visualizer.fit(X_train, y_train)        
	visualizer.score(X_test, y_test)
	visualizer.finalize()
	plt.savefig(os.path.join(visualizer_out_path_sub,save_name)) 
	plt.clf()  

def precision_recall_curve(X_train, y_train, X_test, y_test, RMSD_mode, model,save_name,run_name=None, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs", run_name)):
		os.makedirs(os.path.join(run_dir,"graphs", run_name))
	visualizer_out_path_sub=os.path.join(run_dir,"graphs", run_name)
	if isinstance(y_train, (np.ndarray)) ==False:
		y_train=y_train.to_numpy()
	class_labels=np.unique(y_train).tolist()
	visualizer = PrecisionRecallCurve(model, classes=class_labels,
		colors=["purple", "cyan", "blue"],
		iso_f1_curves=True,
		per_class=True,
		micro=False)
	if type(model).__name__=="XGBoost":
				KX_train, KX_test, Ky_train, Ky_test = KX_train.to_numpy(), KX_test.to_numpy(), Ky_train.to_numpy(), Ky_test.to_numpy()
				num_class=len(np.unique(Ky_train))
				param={'objective':'binary:logistic'}
				if num_class>2:
					Ky_train=Ky_train-1
					Ky_test=Ky_test-1
					param['objective']='multi:softmax'
					param['num_class']=num_class
	visualizer.fit(X_train, y_train)        
	visualizer.score(X_test, y_test)
	visualizer.finalize()
	plt.savefig(os.path.join(visualizer_out_path_sub,save_name)) 
	plt.clf()  
	
def learning_curve(X_train, y_train, model,save_name,  run_name=None, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs", run_name)):
		os.makedirs(os.path.join(run_dir,"graphs", run_name))
	visualizer_out_path_sub=os.path.join(run_dir,"graphs", run_name)

	kfold = KFold(n_splits=10, shuffle=True)
	sizes = np.linspace(0.3, 1.0, 10)
	mcc_scorer=make_scorer(matthews_corrcoef)
	visualizer = LearningCurve(
								model, cv=kfold, scoring=mcc_scorer, train_sizes=sizes, n_jobs=1)

	visualizer.fit(X_train, y_train)
	visualizer.finalize()
	plt.savefig(os.path.join(visualizer_out_path_sub,save_name))
	plt.clf()  
####################################################################################################################
#FEATURE SELECTION
##################

def random_forest_features(X_train, y_train, RMSD_mode, save_name,run_name=None, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs", run_name)):
		os.makedirs(os.path.join(run_dir,"graphs", run_name))
	visualizer_out_path_sub=os.path.join(run_dir,"graphs", run_name)

	# fitting the model
	model = RandomForestClassifier(n_estimators=500, random_state=42)
	model.fit(X_train, y_train)# plotting feature importances
	features = X_train.columns
	importances = model.feature_importances_
	indices = np.argsort(importances)
	plt.figure(figsize=(15,20))
	plt.title('Feature Importances')
	plt.barh(range(len(indices)), importances[indices], color='b', align='center')
	plt.yticks(range(len(indices)), [features[i] for i in indices])
	plt.xlabel('Relative Importance')
	plt.savefig(os.path.join(visualizer_out_path_sub,save_name))
	plt.clf()  

def feature_importance(X_train, y_train,model,save_name,run_name=None, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs", run_name)):
		os.makedirs(os.path.join(run_dir,"graphs", run_name))
	visualizer_out_path_sub=os.path.join(run_dir,"graphs", run_name)

	try:
		labels=X_train.columns
		viz = FeatureImportances(model, labels=labels, relative=False) 
		viz.fit(X_train, y_train)
		viz.show(outpath=os.path.join(visualizer_out_path_sub,save_name)) 
		viz = FeatureImportances(model, labels=labels, stack=True, relative=False)
		viz.fit(X_train, y_train)
		viz.finalize()
		plt.savefig(os.path.join(visualizer_out_path_sub,"stacked_"+save_name)) 
		plt.clf()  
	except:
		try:
			for alg in model.named_estimators:
				model_alg = model.named_estimators[alg]
				viz = FeatureImportances(model_alg, labels=labels, relative=False) 
				viz.fit(X_train, y_train)
				viz.show(outpath=os.path.join(visualizer_out_path_sub,"voting_"+str(alg)+save_name)) 
				viz = FeatureImportances(model_alg, labels=labels, stack=True, relative=False)
				viz.fit(X_train, y_train)
				viz.finalize()
				plt.savefig(os.path.join(visualizer_out_path_sub,"voting_"+str(alg)+"stacked_"+save_name))
				plt.clf()  
		except:
			pass

def recursive_feature_elimination(X_train, y_train,model,save_name,  run_name=None, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs", run_name)):
		os.makedirs(os.path.join(run_dir,"graphs", run_name))
	visualizer_out_path_sub=os.path.join(run_dir,"graphs", run_name)

	kfold = KFold(n_splits=10, shuffle=True)
	mcc_scorer=make_scorer(matthews_corrcoef)
	try:
		X_normalized = preprocessing.normalize(X_train, norm='l2')
		y_train=y_train.astype(str)
		estimator = SVC()
		estimator=estimator.fit(X_train,y=y_train)
		visualizer = RFECV(model, cv=kfold, scoring=mcc_scorer)
		visualizer.fit(X_train, y_train)
		visualizer.finalize()
		plt.savefig(os.path.join(visualizer_out_path_sub,save_name))
		plt.clf()  
	except:
		pass

def xgb_feature_importance():

	xgb.importance(feature_names = colnames(sparse_matrix), model = bst)
####################################################################################################################
#MODEL OPTIMIZATION
###################
def validation_curve(X_train, y_train, RMSD_mode, model,save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
		os.makedirs(os.path.join(run_dir,"graphs"))
	mcc_scorer=make_scorer(matthews_corrcoef)
	param_name=""
	param_range=""
	viz = ValidationCurve(
			model, param_name=param_name,
			param_range=param_range, cv=10, scoring=mcc_scorer)

	viz.fit(X_train, y_train)
	viz.finalize()
	plt.savefig(os.path.join(run_dir,"graphs",save_name)) 
	plt.clf()  


def opti(MCC_df, MCC_df_non_optimized, models_torun):
	"""
	Input: MCC_dic_opti ---  Dictionary containing two key-value pairs per model structure. 
							1) "Model_name":model name
							2) "MCC_mean": mean of all MCC values collected for 10 folds
			MCC_dic  --- same as MCC_dic_opti but with the non-optimized models
			models_torun  --- single-row dataframe with information on model-type, 
								model performance, what the target class is and the model itself.
	"""
	pass

#####################################################################################################################
#MODEL ERRORS
#############
def num_error_gradient_three(metrics,y_test,save_name, run_dir=None):
	if not os.path.exists(os.path.join(run_dir,"graphs")):
            os.makedirs(os.path.join(run_dir,"graphs"))

	percent_in_bounds = metrics['in_bounds'].mean() * 100
	metrics_to_plot = metrics[[c for c in metrics if 'absolute_error' in c]]

	cols=metrics_to_plot.columns
	
	ax = sns.boxplot(data=metrics_to_plot)
	plt.title("Error Plot for Range Numerical Predictor")
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()

	sns.set(font_scale=1)
	fig=sns.distplot(metrics["absolute_error_lower"])
	plt.xlabel(" RMSD in Angstrom")
	plt.ylabel("Nr. of Models")
	plt.title("Error in RMSD in Angstrom - Lower Numerical Predictor")
	plt.savefig(os.path.join(run_dir,"graphs","lower_dist_"+save_name))
	plt.clf()

	sns.set(font_scale=1)
	fig=sns.distplot(metrics["absolute_error_upper"])
	plt.xlabel(" RMSD in Angstrom")
	plt.ylabel("Nr. of Models")
	plt.title("Error in RMSD in Angstrom - Upper Numerical Predictor")
	plt.savefig(os.path.join(run_dir,"graphs","upper_dist_"+save_name))
	plt.clf()

	sns.set(font_scale=1)
	fig=sns.distplot(metrics["absolute_error_mid"])
	plt.xlabel(" RMSD in Angstrom")
	plt.ylabel("Nr. of Models")
	plt.title("Error in RMSD in Angstrom - Mid Numerical Predictor")
	plt.savefig(os.path.join(run_dir,"graphs","mid_dist_"+save_name))
	plt.clf()


	plot=JointPlotVisualizer(size=(1620, 1080), correlation="pearson")
	plot.fit_transform(X=y_test,y=metrics[["absolute_error_mid"]].to_numpy().flatten())
	plot.finalize()
	plt.ylabel("RMSD in Angstroms")
	plt.xlabel("RMSD Error in Angstroms")
	plt.title("Error in RMSD in Angstrom vs RMSD")
	plt.savefig(os.path.join(run_dir,"graphs","Error_vs_RMSD"+save_name))
	plt.clf()




if __name__ == '__main__':
	ROCAUC_curve(None,None,None,None,None,None,None)