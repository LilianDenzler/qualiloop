#!/usr/bin/env python
import os
import numpy as np
import pandas as pd



def make_logo(save_dir):
	df=pd.DataFrame(columns=["pos","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"])
	if not os.path.isdir(os.path.join(save_dir,"aa_dist_CDRH3")):
					os.makedirs(os.path.join(save_dir,"aa_dist_CDRH3"))
	for file in os.listdir(os.path.join(save_dir,"aa_dist_CDRH3")):
		filepath=os.path.join(os.path.join(save_dir,"aa_dist_CDRH3",file))
		pos=int(file.replace(".csv", ""))
		with open(filepath) as file:
			dic={}
			for line in file:
				if not line.startswith("#"):
					line=line.split(",")
					print(line)
					aa=line[0]
					if aa=="-":
						continue
					weight=line[2].replace("\n", "")
					dic[aa]=weight
				dic["pos"]=pos
				
			
			df=df.append(dic, ignore_index=True)
				
	df=df.set_index('pos')
	df=df.sort_index()
	print(df)
	df.to_csv(os.path.join(save_dir,"aa_dist_CDRH3","pssm_CDRH3.csv"), sep='	')