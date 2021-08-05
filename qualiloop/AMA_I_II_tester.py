#!/usr/bin/python3
"""/*************************************************************************
   Program:    Antibody Modelling Assessment tester
   File:       AMA_I_II_tester.py
   
   Version:    V1.1
   Date:       30.04.21
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
   This program takes the models published in the AMA-I and AMA-II studies 
   conducted by Almagro et. al. The RMSD values of the CDRH3 loop are calculated
   for all of the published modelling methods, as well as for abYmod. Re-calculating 
   the CDR-H3 RMSD values for all ensures that the values are comparable. This
   is importatn as teh primary aim of the program is to determine how well 
   abYmod modells the CDR-H3 loop compared to the other programs. 

**************************************************************************
   Usage:
   ======
   This library is intended for analysing abYmod
   . 
**************************************************************************
   Revision History:
   =================
   V0.1   27.04.21 Original
*************************************************************************/"""
import sys
import os
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./')
sys.path.append('~/sync_project/WWW/CDRH3loop')
import subprocess
import pandas as pd
import numpy as np
import joblib
import seaborn as sns
import shutil
import modelmaker
import matplotlib.pyplot as plt
import glob


def get_file_paths():
	pdbabnum_path=glob.glob('**/libs/**/pdbabnum', recursive=True)
	if len(pdbabnum_path)>1:
		raise ValueError('more than one pdbabnum file was found')
	elif len(pdbabnum_path)==0:
		pass
		#raise ValueError('no pdbabnum file was found')
	pdbabnum_path="pdbabnum"

	pdbhstrip_path=glob.glob('**/libs/**/pdbhstrip', recursive=True)
	if len(pdbhstrip_path)>1:
		raise ValueError('more than one pdbhstrip file was found')
	elif len(pdbhstrip_path)==0:
		pass
		#raise ValueError('no pdbhstrip file was found')
	pdbhstrip_path="pdbhstrip"
	pdb2pir_path=glob.glob('**/libs/**/pdb2pir', recursive=True)
	if len(pdb2pir_path)>1:
		print(pdbsolv_path)
		raise Warning('more than one pdb2pir file was found')
	elif len(pdb2pir_path)==0:
		pass
		#raise ValueError('no pdb2pir file was found')
	pdb2pir_path="pdb2pir"

	return pdbabnum_path,pdbhstrip_path,pdb2pir_path

def actual_pirs(actual_dir,pdbabnum_path,pdbhstrip_path):
	for actual_file in os.listdir(actual_dir):
		pdb_ID=str(actual_file.replace(".pdb",""))

		actual_file=os.path.join(actual_dir, actual_file)
		command="{} -c {}".format(pdbabnum_path,actual_file)
		pdbabnum_out=subprocess.check_output(command, shell=True)
		with open(actual_file, "w") as f:
			pdbabnum_out=str(pdbabnum_out).split("\\n")
			pdbabnum_out=["%s\n" % l for l in pdbabnum_out]
			for i in pdbabnum_out:
				if "ATOM" in i or "CONNECT" in i or "MASTER" in i or "END" in i:
					f.write(i.replace("', '", ""))
				
		command="{} {}".format(pdbhstrip_path,actual_file)
		pdbhstrip_out=subprocess.check_output(command, shell=True)
		with open(actual_file, "w") as f:
			pdbhstrip_out=str(pdbabnum_out).split("\\n")
			pdbhstrip_out=["%s\n" % l for l in pdbhstrip_out]
			for i in pdbhstrip_out:
				if "ATOM" in i or "CONNECT" in i or "MASTER" in i or "END" in i:
					f.write(i.replace("', '", ""))

def check_seq(file, actual_file,pdb2pir_path):
	command="{} {} {}".format(pdb2pir_path,file, file+".pir")
	pdb2pir_out_file=subprocess.call(command, shell=True)
	with open(file+".pir","r") as f:
		file_pir=f.read()

	command="{} {} {}".format(pdb2pir_path,actual_file, actual_file+".pir")
	pdb2pir_out_actual=subprocess.call(command, shell=True)
	with open(actual_file+".pir","r") as f:
		actual_file_pir=f.read()
	if actual_file_pir==file_pir:
		pass
	else:
		print(actual_file_pir)
		print(file_pir)
	
def get_RMSD_vals(model_dir, actual_dir, AMA, method, AMA_dir,pdbabnum_path,pdbhstrip_path,pdb2pir_path):
	full_RMSD_df=pd.DataFrame(columns=["local_AA", "local_CA", "global_AA", "global_CA","ID", "method"])
	for file in os.listdir(model_dir):
		file=os.path.join(model_dir, file)
		
		command="{} -c {}".format(pdbabnum_path,file)
		pdbabnum_out=subprocess.check_output(command, shell=True)
		#print(str(pdbabnum_out, 'utf-8'))
		
		with open(file, "w") as f:
			pdbabnum_out=str(pdbabnum_out).split("\\n")
			pdbabnum_out=["%s\n" % l for l in pdbabnum_out]
			for i in pdbabnum_out:
				if "ATOM" in i or "CONNECT" in i or "MASTER" in i or "END" in i:
					f.write(i.replace("', '", ""))
					
		command="{} {}".format(pdbhstrip_path,file)
		pdbhstrip_out=subprocess.check_output(command, shell=True)
		#print(str(pdbhstrip_out, 'utf-8'))
		with open(file, "w") as f:
			pdbhstrip_out=str(pdbabnum_out).split("\\n")
			pdbhstrip_out=["%s\n" % l for l in pdbhstrip_out]
			for i in pdbhstrip_out:
				if "ATOM" in i or "CONNECT" in i or "MASTER" in i or "END" in i:
					f.write(i.replace("', '", ""))				
				

	for file in os.listdir(model_dir):
		skip_file=False
		#while skip_file==False:
		if ".pir" in file:
			continue
		ID=None
		filename=file.split(".")[0]
		for i in ["C705","X836","101F","B21M","C1068","C706","2507","323B3","169B3","Ab02","Ab03","Ab04","Ab05","Ab06","Ab07","Ab08","Ab09","Ab010","Ab011"]:
			if i in filename:
				ID=i
				break
			if "Ab010" not in filename and "Ab011" not in filename and "Ab01" in filename:
				ID="Ab01"
				break
		if ID==None:
			pass
			#input(filename)

			

		for actual_file in os.listdir(actual_dir):
			actual_filename=actual_file.replace(".pdb", "")
			if actual_filename in file:
				actual_file=os.path.join(actual_dir, actual_file)
				file=os.path.join(model_dir, file)
				print(file, actual_file)
				check_seq(file, actual_file,pdb2pir_path)

				RMSD_num_df=modelmaker.get_RMSD_num(file, actual_file, fast_mode=None)
				ID_df=pd.DataFrame()
				ID_df["ID"]=[ID]
				ID_df["method"]=[method]
				RMSD_num_df = pd.concat([RMSD_num_df, ID_df], axis=1)
				print(RMSD_num_df)

				try:
					test=RMSD_num_df[["global_CA"]]
				except:
					full_RMSD_df = pd.concat([RMSD_num_df, full_RMSD_df], axis=0)
					skip_file=True
					break
				if full_RMSD_df.loc[full_RMSD_df['ID'] == ID].empty==False:
					same_code=full_RMSD_df.loc[full_RMSD_df['ID'] == ID]
					same_val=same_code["global_CA"].to_numpy()[0]
					if float(same_val) > float(RMSD_num_df[["global_CA"]].to_numpy()[0]):
						skip_file=True
						break
					else:
						same_code=RMSD_num_df
						skip_file=True
						break
				else:
					full_RMSD_df = pd.concat([RMSD_num_df, full_RMSD_df], axis=0)
					skip_file=True
					break
						

	full_RMSD_df.to_csv(os.path.join(AMA_dir,"RMSD_vals_"+str(AMA)+".csv"))
	#input(full_RMSD_df)
	visualize_per_model(full_RMSD_df,AMA_dir)
	return full_RMSD_df

def visualize_per_model(full_RMSD_df,AMA_dir):
	global_CA = full_RMSD_df["global_CA"].to_list()
	global_CA=[float(x) for x in global_CA]
	full_RMSD_df["global_CA"]=global_CA
	ax = sns.barplot(x="ID", y="global_CA", data=full_RMSD_df, ci=None)
	ax.set_title('RMSD values of the CDR-H3 region for {} - AMA-I'.format(method))
	ax.set(xlabel="PDB code", ylabel = 'RMSD (Å)')
	ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
	ax.legend(bbox_to_anchor=(5.0, 1.0), loc='upper left')
	fig=ax.get_figure()
	plt.tight_layout()
	plt.savefig(os.path.join(AMA_dir,"bar_chart_2"+str(AMA)+str(method)+".png"))
	plt.clf()
	if method=="abYmod":
		local_CA = full_RMSD_df["local_CA"].to_list()
		local_CA=[float(x) for x in local_CA]
		full_RMSD_df["local_CA"]=local_CA
		ax = sns.barplot(x="ID", y="local_CA", data=full_RMSD_df, ci=None)
		ax.set_title('RMSD values of the CDR-H3 region for {} - AMA-I'.format(method))
		ax.set(xlabel="PDB code", ylabel = 'Local Calpha RMSD (Å)')
		ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
		ax.set(xlabel="PDB code", ylabel = 'RMSD (Å)')
		fig=ax.get_figure()
		plt.tight_layout()
		plt.savefig(os.path.join(AMA_dir,"bar_chart_2"+str(AMA)+str(method)+"local_CA"+".png"))
		plt.clf()


def visualize_all(full_RMSD_df, AMA_dir):
	global_CA = full_RMSD_df["global_CA"].to_list()
	global_CA=[float(x) for x in global_CA]
	full_RMSD_df["global_CA"]=global_CA
	ax = sns.barplot(x="ID", y="global_CA", data=full_RMSD_df, hue="method",ci=None)
	ax.set_title('RMSD values of the CDR-H3 region - AMA-I')
	ax.set(xlabel="PDB code", ylabel = 'RMSD (Å)')
	ax.set_xticklabels(ax.get_xticklabels(),rotation=45)
	ax.set(xlabel="PDB code", ylabel = 'RMSD (Å)')
	fig=ax.get_figure()
	plt.tight_layout()
	plt.savefig(os.path.join(AMA_dir,"bar_chart_2"+str(AMA)+"ALL"+".png"))
	plt.clf()
if __name__ == '__main__':
	#python3.6 myfunctions.py ~/sync_project/AMA_I/AMAI_abYmod/ ~/sync_project/AMA_I/actual_PDBs_AMA_I/ ~/sync_project/Redundant_PDBs.txt
	pdbabnum_path,pdbhstrip_path,pdb2pir_path=get_file_paths()
	AMA_I_dir=sys.argv[1]
	actual_dir=os.path.join(AMA_I_dir,"actual_PDBs_AMA_I")
	
	AMA="AMA_I"
	actual_pirs(actual_dir,pdbabnum_path,pdbhstrip_path)
	full_RMSD_df=pd.DataFrame()
	for method in ["ACC", "CCG", "PIGS", "ROS","abYmod"]:
		model_dir=os.path.join(AMA_I_dir,"AMAI_"+str(method))
		single_df=get_RMSD_vals(model_dir, actual_dir, AMA, method,AMA_I_dir,pdbabnum_path,pdbhstrip_path,pdb2pir_path)
		full_RMSD_df = pd.concat([single_df, full_RMSD_df], axis=0)
	visualize_all(full_RMSD_df,AMA_I_dir)
	
	AMA_II_dir=sys.argv[2]
	actual_dir=os.path.join(AMA_II_dir,"actual_PDBs_AMA_II")
	AMA="AMA_II"
	actual_pirs(actual_dir,pdbabnum_path,pdbhstrip_path)
	full_RMSD_df=pd.DataFrame()
	for method in ["ACC","CCG", "JEF", "JOA", "MMT", "PIGS", "SCH","abYmod"]:
		model_dir=os.path.join(AMA_II_dir,"AMAII_"+str(method))
		single_df=get_RMSD_vals(model_dir, actual_dir, AMA, method,AMA_II_dir,pdbabnum_path,pdbhstrip_path,pdb2pir_path)
		full_RMSD_df = pd.concat([single_df, full_RMSD_df], axis=0)
	visualize_all(full_RMSD_df,AMA_II_dir)
	

#note: for the AMAII structure Ab010.pdb the chain A was taken as L-chain and the chain B was taken as H-chain for abymod
#for AMAII structure Ab01.pdb the RMSD could not be calculated for models of the programs SCH, PIGS and MMT.

#note for C705 in AMAI, chains A and B were removed for abymod modelling, H and aL were kept and used
