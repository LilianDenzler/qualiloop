#!/usr/bin/python3

#########################################################################################
#Purpose of Programme:
#1: take each modelled loop structure and run loop flapper programme
#2: calculate different features of each flapped structure: i.e. ecalc (only VdW, electrostats,..), accessibility, packing quality, etc. 
#   build machine learning model 
#3: output the optimal angle of the structure
#4: assess quality of loop flapper programm by calculating the angle of the original PDB file and comparing with prediction of optimal angle


import sys
import os
import signal
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./CDRH3lib')

import numpy as np
import csv
import pandas as pd
import itertools
import math
from biopandas.pdb import PandasPdb
from scipy.spatial.transform import Rotation
import subprocess
import qualiloop
from qualiloop import myfunctions
from functools import reduce
import shutil
#PART 1: Reading the input PDB model

def input_process(pdb_file):
	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb=ppdb.df['ATOM']

	#get index of start-atom and stop-atom of loop
	index_95 = ppdb.loc[(ppdb['residue_number'] == 95)& (ppdb['chain_id']=="H")].index.tolist()[0]
	index_102 = ppdb.loc[(ppdb['residue_number'] == 102)&(ppdb['chain_id']=="H")].index.tolist()[-1]


	#get coordinates of full structure
	full_structure=ppdb[["x_coord", "y_coord", "z_coord"]]
	full_structure=full_structure.to_numpy()

	#get base 95 point
	ppdb95=ppdb.loc[ppdb['residue_number'] == 95]
	ppdb95=ppdb95.loc[ppdb['chain_id'] == "H"]
	ppdb95=ppdb95.loc[ppdb['atom_name'] == "CA"]
	ppdb95=ppdb95[["x_coord","y_coord", "z_coord"]]
	ppdb95=ppdb95.to_numpy()

	#get base 102
	ppdb102=ppdb.loc[ppdb['residue_number'] == 102]
	ppdb102=ppdb102.loc[ppdb['chain_id'] == "H"]
	ppdb102=ppdb102.loc[ppdb['atom_name'] == "CA"]
	ppdb102=ppdb102[["x_coord","y_coord", "z_coord"]]
	ppdb102=ppdb102.to_numpy()
	return (ppdb, index_95,index_102, full_structure, ppdb95, ppdb102)



#shift the whole structure so that the base point H95 is at the origin; then rotate so second basepoint H102 as on x-axis.
def move_to_origin(full_structure, ppdb95, ppdb102):
	#move entire structure so that resH95 is at origin
	for i in range(0, np.shape(full_structure)[0]):
		#for x
		full_structure[i, 0]=full_structure[i, 0]-ppdb95[0,0]
		#for y
		full_structure[i, 1]=full_structure[i, 1]-ppdb95[0,1]
		#for z
		full_structure[i, 2]=full_structure[i, 2]-ppdb95[0,2]
	#define new H102 base-point position after transformation
	ppdb102[0,0]=ppdb102[0,0]-ppdb95[0,0]
	ppdb102[0,1]=ppdb102[0,1]-ppdb95[0,1]
	ppdb102[0,2]=ppdb102[0,2]-ppdb95[0,2]

	#align vectors from origin to res H102 and vector from origin along x axis
	rotation_origin, rmsd=Rotation.align_vectors(ppdb102,[[1,0,0]])
	rotated_full_structure=rotation_origin.apply(full_structure)
	rotated_full_df = pd.DataFrame(data=rotated_full_structure, columns=["x_coord", "y_coord", "z_coord"])

	return(ppdb95, ppdb102, rotated_full_df)


#PART2 flap the loop and make PDB files with altered angles

#rotate the loop structure by one degree around the x-axis
#outputs the rotated_full_df that has a changed loop structure, flapped by 1 degree
def flap_loop(rotated_full_df, index_95,index_102, angle_diff):
	r = Rotation.from_euler('x', angle_diff, degrees=True)
	for i in range(index_95, index_102+1):
		a_loop_vector=rotated_full_df.iloc[[i]]
		new_loop_vector=r.apply(a_loop_vector)[0]
		rotated_full_df.loc[i] = new_loop_vector

	return (rotated_full_df)

#Creates a PDB file in the save_path directory; name of file will be flapped_structure and the difference in angle to the original appended. 
#takes the structure at it's original position with the rotated loop as input
def pdb_integrate(original_rotation_structure, pdb_file, angle_diff, save_path):
	new_structure=pd.DataFrame(data=original_rotation_structure, index=None, columns=["x_coord","y_coord","z_coord"])
	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb.df['ATOM'][['x_coord','y_coord','z_coord']] = new_structure
	path=os.path.join(save_path,"flapped_structure_"+str(angle_diff)+".pdb")
	ppdb.to_pdb(path=path, 
            records=None, 
            gz=False, 
            append_newline=True)
	return (None)

#full loop flapping process, flap the loop from 
def full_flap (rotated_full_df, index_95, index_102, rotation_structure, ppdb95, pdb_file, save_path, max_angle_diff):
	for direction in ["clockwise", "counterclockwise"]:
		if direction=="clockwise":
			angles=range(0, int(max_angle_diff)+1,1)
		else:
			angles=range(0,(-1)*int(max_angle_diff)-1,-1)
		for angle_diff in angles:
			flapped_full_df=flap_loop(rotated_full_df, index_95,index_102, angle_diff)
			#return to initial position
			#first reverse rotation to place H102 on x-axis
			reverse_rotation=rotation_structure.inv()
			original_rotation_structure=reverse_rotation.apply(flapped_full_df)
			#then reverse transformation to place H95 at origin
			for i in range(0, np.shape(original_rotation_structure)[0]):
				#for x
				original_rotation_structure[i, 0]=original_rotation_structure[i, 0]+ppdb95[0,0]
				#for y
				original_rotation_structure[i, 1]=original_rotation_structure[i, 1]+ppdb95[0,1]
				#for z
				original_rotation_structure[i, 2]=original_rotation_structure[i, 2]+ppdb95[0,2]
			pdb_integrate(original_rotation_structure, pdb_file, angle_diff, save_path)
	return None

def run(pdb_file, save_path, max_angle_diff):
	#process the input pdb model file, get the index of the base residues
	(ppdb, index_95,index_102, full_structure, ppdb95, ppdb102)=input_process(pdb_file)

	#move the full structure so that H95 is on the origin and H102 is on the x-axis 
	(ppdb95, ppdb102, rotated_full_df)=move_to_origin(full_structure, ppdb95, ppdb102)
	
	#align vectors from origin to res H102 and vector from origin along x axis
	rotation_structure, rmsd=Rotation.align_vectors(ppdb102,[[1,0,0]])
	rotated_full_structure=rotation_structure.apply(full_structure)
	rotated_full_df = pd.DataFrame(data=rotated_full_structure, columns=["x_coord", "y_coord", "z_coord"])

	#flap the loop, output pdb files for all angles from -max_angle_diff to +max_angle_diff; 
	#file of 0 degree angle difference should correspond to original file, use as sanity check
	full_flap (rotated_full_df, index_95, index_102, rotation_structure, ppdb95, pdb_file, save_path, max_angle_diff)
	return None

def get_ecalc(save_path):
	os.environ['ECALCDATA'] = '/serv/www/html_lilian/libs/ecalc-master/data'
	os.environ['BIOPLIB_DEPRECATED_QUIET']='silence'
	os.environ['DATADIR']='/serv/www/html_lilian/libs/data'
	signal.signal(signal.SIGSEGV, signal.SIG_IGN)
	energy_list=[]
	for file in os.listdir(save_path):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(save_path,file)
			command='pdbgetchain L,H {}|/serv/www/html_lilian/libs/mutmodel -R | /serv/www/html_lilian/libs/bioptools-master/src/pdbchain | /serv/www/html_lilian/libs/bioptools-master/src/pdbhadd -c | /serv/www/html_lilian/libs/bioptools-master/src/pdbcter -c | /serv/www/html_lilian/libs/bioptools-master/src/pdbrenum > {}'.format(file,os.path.join(save_path, filename+'.pdh'))
			p = subprocess.call(command, shell=True)
			command='/serv/www/html_lilian/libs/ecalc-master/src/ecalc -p {}'.format(os.path.join(save_path, filename+'.pdh'))
			ecalc_out=subprocess.check_output(command, shell=True)
			ecalc_out=str(ecalc_out, 'utf-8')	
					
			for line in ecalc_out.split('\n'):
				if 'Base structure energy' in line:
					energy=line.split(' = ')[1]
					energy=energy.replace('/n','')
					energy=float(energy)
					energy_list+=[[str(file), energy]]
	
	energy_df = pd.DataFrame(data=energy_list, columns=['file', 'energy'], index=None)
	return energy_df

def accessibility (save_path):
	accessibility_df=pd.DataFrame()
	for file in os.listdir(save_path):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(save_path,file)
			if accessibility_df.empty != True:
				new_line=myfunctions.accessibility(file)
				new_line["file"]=str(file)
				accessibility_df=pd.concat([new_line, accessibility_df])
			else:
				accessibility_df=myfunctions.accessibility(file)
				accessibility_df["file"]=str(file)
	return (accessibility_df)

def get_protrusion(pdb_file):
	protrusion_df= myfunctions.get_protrusion(pdb_file)
	protrusion_res= protrusion_df["tip"].loc[0][1:]
	protrusion_res=int(protrusion_res)

	ppdb = PandasPdb().read_pdb(os.path.join(pdb_file))
	ppdb=ppdb.df['ATOM']
	protrusion_res = ppdb.loc[(ppdb['residue_number'] == protrusion_res)& (ppdb['chain_id']=="H")& (ppdb['atom_name'] == "CA")]
	A=protrusion_res[["x_coord","y_coord","z_coord"]]
	A=A.to_numpy()[0]
	return A

def get_original_angle(protrusion_point):
	#A=maximal protruion point, point B is defined as point A projected onto the Y-Z Plane
	#As the rotation is only about the X-axis, the X-coordinate can be disregarded
	B=[0,0,0]
	angle=math.atan2(protrusion_point[2], protrusion_point[1]) - math.atan2(B[2], B[1])
	return angle

def get_angle_df(save_path):
	angle_df=pd.DataFrame()
	for file in os.listdir(save_path):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(save_path,file)
			if angle_df.empty != True:
				A=get_protrusion(file)
				angle=get_original_angle(A)
				new_line=pd.DataFrame(data={"angle":[angle]})
				new_line["file"]=str(file)
				angle_df=pd.concat([new_line, angle_df])
			else:
				A=get_protrusion(file)
				angle=get_original_angle(A)
				angle_df=pd.DataFrame(data={"angle":[angle]})
				angle_df["file"]=str(file)
				
	return (angle_df)



def merge_features(save_path):
	accessibility_df=accessibility (save_path)
	print("ACCESSIBILITY;", accessibility_df)
	energy_df=get_ecalc(save_path)
	print("ENERGY;",energy_df)
	angle_df=get_angle_df(save_path)
	print("ANGLE;", angle_df)
	if accessibility_df.empty!=False or energy_df.empty!=False or angle_df.empty!=False:
		df_final=pd.DataFrame()
		df_file_angles=pd.DataFrame()
		return df_final, df_file_angles
	feature_df_list=[accessibility_df, energy_df, angle_df]
	df_final = reduce(lambda left,right: pd.merge(left,right,on='file'), feature_df_list)
	df_file_angles=pd.DataFrame()
	selected_columns=df_final[["file"]]
	df_file_angles = selected_columns.copy()
	df_file_angles["angle"]=df_final[["angle"]]
	df_final.drop(labels="file",axis=1)
	return df_final, df_file_angles

def df_to_row(df_final, max_angle_diff):
	column_names=[]
	column_row=[]
	for i in df_final.columns:
		for a in range(((-1)*int(max_angle_diff)), (int(max_angle_diff)+1)):
			column_names+=[i+str(a)]
		add_columns=df_final[[i]].to_numpy()
		np.append(column_row, add_columns)
	row_df=pd.DataFrame(data=column_row, columns=column_names, index=None)
	return row_df

#do ML prediction using this feature table
def maḱe_train_test_df(actual_pdbs, save_path, max_angle_diff):
	full_feature_df=pd.DataFrame()
	for file in os.listdir(actual_pdbs):
		if file.endswith(".pdb"):
			filename=os.path.splitext(os.path.basename(file))[0]
			file=os.path.join(actual_pdbs,file)
		else:
			continue

		print("FILENAME:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!",filename)
		#create loop flapping structures for file
		if not os.path.exists(save_path):
			os.mkdir(save_path)
		run(file, save_path, max_angle_diff)
		#make individual df of this file 
		df_final, df_file_angles=merge_features(save_path)
		if df_final.empty!=False or df_file_angles.empty!=False:
			continue
		#compress df into single row df
		row_df=df_to_row(df_final, max_angle_diff)
		if full_feature_df.empty != True:
			full_feature_df=pd.concat([full_feature_df, row_df])
		else:
			full_feature_df=row_df
	shutil.rmtree(save_path)
	msk = np.random.rand(len(full_feature_df)) < 0.8
	train_df = full_feature_df[msk]
	test_df = full_feature_df[~msk]

	return train_df, test_df


def train_model(train_df, test_df):
	df_actual=train_df[["angle"]].to_numpy()
	train_y=df_actual
	df_vals=train_df.drop(labels="angle",axis=1).to_numpy()
	train_X=df_vals
	model = RandomForestClassifier(n_estimators=200)
	model.fit(train_X, train_y)
	RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
		max_depth=None, max_features='auto', max_leaf_nodes=None,
		min_impurity_decrease=0.0, min_impurity_split=None,
		min_samples_leaf=1, min_samples_split=2,
		min_weight_fraction_leaf=0.0, n_estimators=200, n_jobs=None,
		oob_score=False, random_state=None, verbose=0,
		warm_start=False)

	num_folds = 10
	acc_per_fold = []
	MCC_per_fold=[]
	fold_no = 1
	kfold = KFold(n_splits=num_folds, shuffle=True)
	fold_dic={}
	for train, test in kfold.split(train_X, train_y):
		print('------------------------------------------------------------------------')
		print(f'Training for fold {fold_no} ...')
		model = RandomForestClassifier(n_estimators=200)
		model.fit(train_X[train], train_y[train])
		RandomForestClassifier(bootstrap=True, class_weight=None, criterion='gini',
			max_depth=None, max_features='auto', max_leaf_nodes=None,
			min_impurity_decrease=0.0, min_impurity_split=None,
			min_samples_leaf=1, min_samples_split=2,
			min_weight_fraction_leaf=0.0, n_estimators=200, n_jobs=None,
			oob_score=False, random_state=None, verbose=0,
			warm_start=False)

		rf_prediction = model.predict(train_X[test])
		MCC_of_fold=matthews_corrcoef(train_y[test],rf_prediction, sample_weight=None )
		scores = model.score(train_X[test], train_y[test])
		acc_per_fold.append(scores * 100)
		MCC_per_fold.append(MCC_of_fold)
		fold_dic.update({MCC_of_fold:fold_no})
		fold_no = fold_no + 1
	# == Provide average scores ==
	print('Score per fold')
	for i in range(0, len(acc_per_fold)):
	  print(f'> Fold {i+1} -  Accuracy: {acc_per_fold[i]}%')
	print('Average scores for all folds:')
	print(f'> Accuracy: {np.mean(acc_per_fold)} (+- {np.std(acc_per_fold)})')
	return (np.mean(MCC_per_fold),rf_prediction)




if __name__ == '__main__':
	actual_pdbs=sys.argv[1]
	save_path=sys.argv[2]
	max_angle_diff=sys.argv[3]

	train_df, test_df=maḱe_train_test_df(actual_pdbs, save_path, max_angle_diff)
	train_model(train_df, test_df)






"""FILENAME:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2WUC_1
ACCESSIBILITY;     Access   Relacc    Scacc  Screlacc  Access_avg  Relacc_avg  Scacc_avg  Screlacc_avg                                           file
0  304.985  152.715  272.533   164.202   38.123125   19.089375  34.066625     20.525250   ./flapped_full_test/flapped_structure_-4.pdb
0  200.181  100.423  179.056   107.426   25.022625   12.552875  22.382000     13.428250    ./flapped_full_test/flapped_structure_1.pdb
0  311.697  151.379  276.544   158.413   38.962125   18.922375  34.568000     19.801625   ./flapped_full_test/flapped_structure_-1.pdb
0  205.414  102.948  184.405   110.045   25.676750   12.868500  23.050625     13.755625    ./flapped_full_test/flapped_structure_2.pdb
0  198.990   99.561  177.582   106.084   24.873750   12.445125  22.197750     13.260500  ./flapped_full_test/flapped_structure_-10.pdb
0  307.757  152.378  274.596   162.238   38.469625   19.047250  34.324500     20.279750   ./flapped_full_test/flapped_structure_-3.pdb
0  260.552  132.455  234.320   145.275   32.569000   16.556875  29.290000     18.159375   ./flapped_full_test/flapped_structure_-7.pdb
0  224.186  112.187  202.630   119.764   28.023250   14.023375  25.328750     14.970500    ./flapped_full_test/flapped_structure_4.pdb
0  312.420  151.203  276.961   157.914   39.052500   18.900375  34.620125     19.739250   ./flapped_full_test/flapped_structure_10.pdb
0  291.503  149.291  259.776   161.041   36.437875   18.661375  32.472000     20.130125   ./flapped_full_test/flapped_structure_-5.pdb
0  272.092  141.563  242.472   154.500   34.011500   17.695375  30.309000     19.312500   ./flapped_full_test/flapped_structure_-6.pdb
0  261.671  133.535  235.147   146.584   32.708875   16.691875  29.393375     18.323000    ./flapped_full_test/flapped_structure_7.pdb
0  280.390  145.607  249.428   158.002   35.048750   18.200875  31.178500     19.750250    ./flapped_full_test/flapped_structure_8.pdb
0  252.074  125.013  229.903   138.007   31.509250   15.626625  28.737875     17.250875    ./flapped_full_test/flapped_structure_6.pdb
0  213.274  106.501  192.799   113.938   26.659250   13.312625  24.099875     14.242250    ./flapped_full_test/flapped_structure_3.pdb
0  312.420  151.203  276.961   157.914   39.052500   18.900375  34.620125     19.739250    ./flapped_full_test/flapped_structure_0.pdb
0  236.716  118.067  213.425   126.871   29.589500   14.758375  26.678125     15.858875    ./flapped_full_test/flapped_structure_5.pdb
0  224.186  112.187  202.630   119.764   28.023250   14.023375  25.328750     14.970500   ./flapped_full_test/flapped_structure_-9.pdb
0  304.985  152.715  272.533   164.202   38.123125   19.089375  34.066625     20.525250    ./flapped_full_test/flapped_structure_9.pdb
0  311.211  152.594  276.498   160.619   38.901375   19.074250  34.562250     20.077375   ./flapped_full_test/flapped_structure_-2.pdb
0  249.548  123.984  226.208   134.901   31.193500   15.498000  28.276000     16.862625   ./flapped_full_test/flapped_structure_-8.pdb
Error: (mutmodel) 
790 hydrogens were added.
Traceback (most recent call last):
  File "flapper_full.py", line 355, in <module>
    train_df, test_df=maḱe_train_test_df(actual_pdbs, save_path, max_angle_diff)
  File "flapper_full.py", line 281, in maḱe_train_test_df
    df_final, df_file_angles=merge_features(save_path)
  File "flapper_full.py", line 237, in merge_features
    energy_df=get_ecalc(save_path)
  File "flapper_full.py", line 165, in get_ecalc
    ecalc_out=subprocess.check_output(command, shell=True)
  File "/usr/lib64/python3.6/subprocess.py", line 356, in check_output
    **kwargs).stdout
  File "/usr/lib64/python3.6/subprocess.py", line 438, in run
    output=stdout, stderr=stderr)
subprocess.CalledProcessError: Command './libs/ecalc-master/src/ecalc -p ./flapped_full_test/flapped_structure_-8.pdh' died with <Signals.SIGSEGV: 11>.
[lilian@serv1 html_lilian]$"""



'''5TIH_1
Too many chains to add disulphide topology
Failure in scanning for disulphides
1026 hydrogens were added.
Too many chains to add disulphide topology
Failure in scanning for disulphides'''

'''
6B0A_1
728 hydrogens were added.
Too many chains to add disulphide topology
Failure in scanning for disulphides
731 hydrogens were added.'''