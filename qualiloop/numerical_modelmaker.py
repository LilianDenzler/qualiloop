import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from keras.models import Sequential
from keras.layers import Dense, Activation
from keras.callbacks import ModelCheckpoint
import keras.backend as K
import tensorflow as tf
import os
from scipy import stats
from qualiloop import visualizer
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import matthews_corrcoef
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
import math
from xgboost import XGBRegressor
from sklearn.metrics import accuracy_score
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import RepeatedStratifiedKFold

import joblib


def correct_within_angstrom(absolute_error_df):
	#get percentage of time the prediction is correct within 1 Angstrom/ .5 ANgstrom
	total=absolute_error_df.shape[0]
	within_1=absolute_error_df[absolute_error_df < 1.0].count()[0]
	percentage_below_1=float(within_1/total)	
	within_05=absolute_error_df[absolute_error_df < 0.5].count()[0]
	percentage_below_05=float(within_05/total)	

	return within_1, percentage_below_1, within_05, percentage_below_05
  

def MCC_of_num_predictor(val_df_mode, model,X_train,y_train,X_test,y_test):
	#calculate the MCC value of the numerical predictor
	upper_thresh_nom=val_df_mode["upper_threshold"]
	lower_thresh_nom=val_df_mode["lower_threshold"]
	nom_label=val_df_mode["nom_label"]
	y_test_nom=0
	y_train_nom=0
	
	try:
		y_test=y_test.numpy()
		y_train=y_train.numpy()
		y_test_nom=y_test.copy()
		y_train_nom=y_train.copy()
	except:
		y_test_nom=y_test.copy()
		y_train_nom=y_train.copy()
	prediction_df=pd.DataFrame()
	prediction_df_train=pd.DataFrame()
	prediction_df['pred'] = model.predict(X_test).flatten()
	prediction_df_train['pred_train'] = model.predict(X_train).flatten()
	for upper,lower,label in zip(upper_thresh_nom,lower_thresh_nom, nom_label):
		prediction_df[(prediction_df['pred']<0)]=None
		prediction_df[(prediction_df['pred'] <=upper) & (prediction_df['pred']>lower)]=label
		prediction_df_train[(prediction_df_train['pred_train'] <=upper) & (prediction_df_train['pred_train']>lower)]=label
		y_train_nom[(y_train_nom <=upper) & (y_train_nom>lower)]=label
		y_test_nom[(y_test_nom <=upper) & (y_test_nom>lower)]=label
	try:
		y_test_nom=y_test_nom.to_numpy()
		y_train_nom=y_train_nom.to_numpy()
	except:
		try:
			y_test_nom=y_test_nom.numpy()
			y_train_nom=y_train_nom.numpy()
		except:

			pass
	if np.any(np.isnan(prediction_df['pred'].to_numpy()))==True:
		print("Predicor predicted values below 0.")

		return None,None,None,None

	prediction_df["y_test"]=y_test_nom
	prediction_df_train["y_train"]=y_train_nom
	wrong_count_test=prediction_df[(prediction_df['pred']!=prediction_df["y_test"])].count()[0]
	wrong_count_train=prediction_df_train[(prediction_df_train['pred_train']!=prediction_df_train["y_train"])].count()[0]
	MCC=matthews_corrcoef(prediction_df["y_test"].to_numpy(),prediction_df['pred'].to_numpy(), sample_weight=None)
	MCC_train=matthews_corrcoef(prediction_df_train["y_train"].to_numpy(),prediction_df_train['pred_train'].to_numpy(), sample_weight=None)
	accuracy = accuracy_score(prediction_df_train["y_train"].to_numpy(),prediction_df_train['pred_train'].to_numpy())
	print("Accuracy: %.2f%%" % (accuracy * 100.0))
	print(prediction_df[['pred',"y_test"]])
	print(y_test_nom)
	print(y_train_nom)


	return wrong_count_test, wrong_count_train, MCC, MCC_train


def gradientBoost_quantile(X_train, y_train, X_test,y_test, val_df_mode,RMSD_name,run_dir=None, already_trained=True):
	if already_trained==False:
		lower_alpha = 0.1
		upper_alpha = 0.9
		lower_model = GradientBoostingRegressor(loss="quantile",                   
											   alpha=lower_alpha)
		mid_model = GradientBoostingRegressor(loss="ls")# ls for mean, GradientBoostingRegressor(loss="quantile", alpha=0.5) for median
		upper_model = GradientBoostingRegressor(loss="quantile",
											   alpha=upper_alpha)

		lower_model.fit(X_train, y_train)
		mid_model.fit(X_train, y_train)
		upper_model.fit(X_train, y_train)
		#save model
		if not os.path.isdir(os.path.join(run_dir,"models")):
			os.makedirs(os.path.join(run_dir,"models"))
		joblib.dump(lower_model, os.path.join(run_dir,"models","lower_model_gradient_boost_regressor.pkl")) 
		joblib.dump(mid_model, os.path.join(run_dir,"models","mid_model_gradient_boost_regressor.pkl")) 
		joblib.dump(upper_model, os.path.join(run_dir,"models","upper_model_gradient_boost_regressor.pkl")) 

	elif already_trained==True:
		if os.path.exists(os.path.join(run_dir,"models","lower_model_gradient_boost_regressor"+".pkl")):
			lower_model=joblib.load(os.path.join(run_dir,"models","lower_model_gradient_boost_regressor"+".pkl"))
		if os.path.exists(os.path.join(run_dir,"models","mid_model_gradient_boost_regressor"+".pkl")):
			mid_model=joblib.load(os.path.join(run_dir,"models","mid_model_gradient_boost_regressor"+".pkl"))
		if os.path.exists(os.path.join(run_dir,"models","upper_model_gradient_boost_regressor"+".pkl")):
			upper_model=joblib.load(os.path.join(run_dir,"models","upper_model_gradient_boost_regressor"+".pkl"))
	predictions = pd.DataFrame()
	predictions["y_test"]=y_test
	predictions['lower'] = lower_model.predict(X_test)
	predictions['mid'] = mid_model.predict(X_test)
	predictions['upper'] = upper_model.predict(X_test)
	predictions['absolute_error_lower'] = (predictions['lower'] - predictions["y_test"]).abs()
	predictions['absolute_error_upper'] = (predictions['upper'] - predictions["y_test"]).abs()
	predictions['absolute_error_interval'] = (predictions['absolute_error_lower'] + predictions['absolute_error_upper']) / 2
	predictions['absolute_error_mid'] = (predictions['mid'] - predictions["y_test"]).abs()
	predictions['in_bounds'] = predictions["y_test"].between(left=predictions['lower'], right=predictions['upper'])
	
	metrics= predictions[['absolute_error_lower', 'absolute_error_upper', 'absolute_error_interval', 'absolute_error_mid', 'in_bounds']].copy()
	visualizer.num_error_gradient_three(metrics,y_test.to_numpy().flatten(),"gradientBoost_quantile.png", run_dir=run_dir)
	within_df=pd.DataFrame(columns=["within_1","percentage_below_1", "within_05", "percentage_below_05"])
	for i in ['absolute_error_lower', 'absolute_error_upper', 'absolute_error_mid']:
		absolute_error_df= predictions[[i]].copy()
		within_1, percentage_below_1, within_05, percentage_below_05=correct_within_angstrom(absolute_error_df)
		within_df.loc[i]=[within_1, percentage_below_1, within_05, percentage_below_05]
	model_MCCs=pd.DataFrame(columns=["wrong_count_test", "wrong_count_train", "MCC", "MCC_train"])
	for model_name,model in zip(["lower_model", "mid_model", "upper_model"],[lower_model, mid_model, upper_model]):
		wrong_count_test, wrong_count_train, MCC, MCC_train=MCC_of_num_predictor(val_df_mode, model, X_train,y_train,X_test,y_test)
		model_MCCs.loc[model_name]=[wrong_count_test, wrong_count_train, MCC, MCC_train]
	sns.jointplot(x=y_test.to_numpy().flatten(),y=mid_model.predict(X_test).flatten(), color=None,kind='reg',stat_func=stats.pearsonr)
	plt.ylabel("Predicted RMSD in Angstroms")
	plt.xlabel("Actual RMSD in Angstroms")
	plt.tight_layout()
	plt.suptitle("Error in RMSD in Angstrom vs RMSD")
	plt.savefig(os.path.join(run_dir,"graphs","gradientBoost_quantile_error.png"))
	plt.clf()

	return predictions, within_df, model_MCCs


def linear_regression_model(X_train, y_train,X_test,y_test,val_df_mode, run_dir=None,already_trained=True):
	if already_trained==False:
		model=LinearRegression()
		model.fit(X_train, y_train)
		if not os.path.isdir(os.path.join(run_dir,"models")):
			os.makedirs(os.path.join(run_dir,"models"))
		joblib.dump(model, os.path.join(run_dir,"models","linear_regression.pkl")) 

	if already_trained==True:
		if os.path.exists(os.path.join(run_dir,"models","linear_regression"+".pkl")):
			model=joblib.load(os.path.join(run_dir,"models","linear_regression"+".pkl"))
			
	linear_regression_preds=pd.DataFrame()
	linear_regression_preds["y_test"]=y_test
	linear_regression_preds["pred"]=model.predict(X_test)
	test_mse = mean_squared_error(y_test, model.predict(X_test))
	test_mae = mean_absolute_error(y_test, model.predict(X_test))
	linear_regression_preds['absolute_error'] = (linear_regression_preds['pred'] - linear_regression_preds["y_test"]).abs()
	linear_regression_description=linear_regression_preds.describe()
	test_rmse=math.sqrt(test_mse)
	#wrong_count_test, wrong_count_train, MCC, MCC_train=MCC_of_num_predictor(val_df_mode, model, X_train,y_train,X_test,y_test)
	return linear_regression_preds,linear_regression_description, test_mse, test_mae, test_rmse#, MCC, MCC_train


def keras_model(X_train, y_train,X_test,y_test,save_name,val_df_mode,run_dir=None,already_trained=True):
	X_train_tf = tf.convert_to_tensor(X_train) 
	y_train_tf = tf.convert_to_tensor(y_train) 
	X_test_tf = tf.convert_to_tensor(X_test) 
	y_test_tf = tf.convert_to_tensor(y_test)
	model = Sequential()
	model.add(Dense(units=3, input_dim=X_train.shape[1],activation='relu'))
	model.add(Dense(units=3, activation='relu'))
	model.add(Dense(1))
	model.compile(loss='mae', optimizer='adadelta')
	print(model.summary())

	if already_trained==True:
		if os.path.exists(os.path.join(run_dir,"models","three_layer"+".ckpt")):
			#model.load_weights(os.path.join(run_dir,"models","three_layer"+".ckpt"))
			keras.models.load_model(os.path.join(run_dir,"models","three_layer"+".ckpt"))


	if already_trained==False:
		model.fit(X_train_tf, y_train_tf, epochs=2000, batch_size=32, verbose=0)
		model.save(os.path.join(run_dir,"models","three_layer"+".ckpt"))
		#cp_callback = ModelCheckpoint(os.path.join(run_dir,"models","three_layer"+".ckpt"),save_weights_only=False, verbose=1)
		#print("saved HERE", os.path.join(run_dir,"models","three_layer"+".ckpt"))


	print(model.evaluate(X_test_tf,y_test_tf))
	y_pred = model.predict(X_test_tf)

	wrong_count_test, wrong_count_train, MCC, MCC_train=MCC_of_num_predictor(val_df_mode, model, X_train,y_train,X_test,y_test)
	sns.jointplot(x=y_test.to_numpy().flatten(),y=y_pred.flatten(), color=None,kind='reg',stat_func=stats.pearsonr)
	plt.ylabel("Predicted RMSD in Angstroms")
	plt.xlabel("Actual RMSD in Angstroms")
	plt.tight_layout()
	plt.suptitle("Error in RMSD in Angstrom vs RMSD")
	plt.savefig(os.path.join(run_dir,"graphs",save_name))
	plt.clf()
	return MCC, MCC_train
	
#https://github.com/sachinruk/KerasQuantileModel/blob/master/Keras%20Quantile%20Model.ipynb
#https://colab.research.google.com/drive/1nXOlrmVHqCHiixqiMF6H8LSciz583_W2#scrollTo=g7s7Grj-A-Sf
############################################################################################
def tilted_loss(q,y,f):
	e = (y-f)
	return K.mean(K.maximum(q*e, (q-1)*e), axis=-1)

def make_model(input_dim):
	model = Sequential()
	model.add(Dense(units=3, input_dim=input_dim,activation='relu'))
	model.add(Dense(units=3, activation='relu'))
	model.add(Dense(1))
	return model

def quantile_loss_keras(X_train, y_train,X_test,y_test, save_name,val_df_mode, run_dir=None, already_trained=False):
	X_train_tf = tf.convert_to_tensor(X_train) 
	y_train_tf = tf.convert_to_tensor(y_train) 
	X_test_tf = tf.convert_to_tensor(X_test) 
	y_test_tf = tf.convert_to_tensor(y_test)
	qs = [0.1, 0.5, 0.9]
	#plt.scatter(X_train,y_train)
	for q,color in zip(qs,["r","g","b"]):
		model = make_model(X_train.shape[1])
		model.compile(loss=lambda y,f: tilted_loss(q,y,f), optimizer='adadelta')
		"""if already_trained==True:
			if os.path.exists(os.path.join(run_dir,"models","three_layer_"+str(q)+".ckpt")):
				keras.models.load_model(run_dir,"models","three_layer_"+str(q)+".ckpt")"""

		if already_trained==False or already_trained==True:
			model.fit(X_train_tf, y_train_tf, epochs=2000, batch_size=32, verbose=0)
			#model.save(os.path.join(run_dir,"models","three_layer_"+str(q)+".ckpt"))
			#don't save model, as unique loss function makes it hard to re-load. Not good model anyways
			
		
		# Predict the quantile
		y_pred = model.predict(X_test_tf)
		print(X_train,y_train,X_test,y_test)
		wrong_count_test, wrong_count_train, MCC, MCC_train=MCC_of_num_predictor(val_df_mode, model, X_train_tf,y_train_tf,X_test_tf,y_test_tf)

		#ax1.scatter(y=y_test, x=y_pred, color=color, label=q) # plot out this quantile
	
		sns.jointplot(y=y_test.to_numpy().ravel(), x=y_pred.ravel(), color=color, label=q,kind='reg',stat_func=stats.pearsonr)
		sns.lmplot(x="predicted", y="actual", data=pd.DataFrame({'predicted':y_pred.ravel(), 'actual':y_test.to_numpy().ravel()}),line_kws={'color': color});
		plt.ylabel("Actual RMSD in Angstroms")
		plt.xlabel("Predicted RMSD in Angstroms")
		plt.suptitle("Predicted vs. Actual RMSD in Angstrom vs RMSD") 
		plt.savefig(os.path.join(run_dir,"graphs",save_name+str(q)+".png"))
		plt.clf()

		return MCC, MCC_train

################################################################################################################

def xgboost_model(X_train, y_train,X_test,y_test, val_df_mode,save_name,run_dir=None, already_trained=False):
	model = XGBRegressor(n_estimators=1000, max_depth=7, eta=0.1, subsample=0.7, colsample_bytree=0.8, reg_lambda=0.8, gamma=24, learning_rate=0.123)
	eval_set = [(X_train, y_train), (X_test, y_test)]

	if already_trained==True:
		if os.path.exists(os.path.join(run_dir,"models","XGBRegressor"+".dat")):
			joblib.load(os.path.join(run_dir,"models","XGBRegressor"+".dat"))
	if already_trained==False:
		model.fit(X_train, y_train, eval_metric=["rmse", "logloss"], eval_set=eval_set, verbose=False)

		joblib.dump(model, os.path.join(run_dir,"models","XGBRegressor"+".dat"))
		
	xgboost_plotter(model,save_name, run_dir=run_dir)
	# make predictions for test data
	y_pred = model.predict(X_test)
	wrong_count_test, wrong_count_train, MCC, MCC_train=MCC_of_num_predictor(val_df_mode, model,X_train,y_train,X_test,y_test)
	return wrong_count_test, wrong_count_train, MCC, MCC_train

def xgboost_plotter(model,save_name, run_dir=None):
	# retrieve performance metrics
	results = model.evals_result()
	epochs = len(results['validation_0']['rmse'])
	x_axis = range(0, epochs)
	# plot log loss
	fig, ax = plt.subplots()
	ax.plot(x_axis, results['validation_0']['logloss'], label='Train')
	ax.plot(x_axis, results['validation_1']['logloss'], label='Test')
	ax.legend()
	plt.ylabel('Log Loss')
	plt.title('XGBoost Log Loss')
	plt.savefig(os.path.join(run_dir,"graphs",save_name+"_log_loss.png"))
	plt.clf()
	# plot classification error
	fig, ax = plt.subplots()
	ax.plot(x_axis, results['validation_0']['rmse'], label='Train')
	ax.plot(x_axis, results['validation_1']['rmse'], label='Test')
	ax.legend()
	plt.ylabel('Root-squared Mean Squared Error')
	plt.title('XGBoost Error')
	plt.savefig(os.path.join(run_dir,"graphs",save_name+"_classification_error.png"))
	plt.clf()


if __name__ == '__main__':
	val_df_mode["upper_threshold"]=[2,4,100]
	val_df_mode["lower_threshold"]=[0,2,4]
	val_df_mode["nom_label"]=[1,2,3]
	print("upper thresholds for categories: ",val_df_mode["upper_threshold"])
	print("lower thresholds for categories: ",val_df_mode["lower_threshold"])
	print("labels for categories: ",val_df_mode["nom_label"])

	xgboost_model(X_train, y_train,X_test,y_test, val_df_mode,"xgboost_regressor",run_dir="./")