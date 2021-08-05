#!/usr/bin/python3
#do chmod +x name_of_script before running

import pandas as pd
import os
import subprocess
import glob

def get_file_paths():
	exposedhphob_path=glob.glob('**/exposedhphob', recursive=True)
	if len(exposedhphob_path)>1:
		raise ValueError('more than one exposedhphob file was found')
	elif len(exposedhphob_path)==0:
		raise ValueError('no exposedhphob file was found')
	bioplib_data_path=glob.glob('**/bioplib-master/data', recursive=True)
	if len(bioplib_data_path)>1:
		raise ValueError('more than one bioplib data file was found')
	elif len(bioplib_data_path)==0:
		raise ValueError('no bioplib data file was found')

	pdbsolv_path=glob.glob('**/pdbsolv', recursive=True)
	if len(pdbsolv_path)>1:
		print(pdbsolv_path)
		#raise Warning('more than one pdbsolv file was found')
	elif len(pdbsolv_path)==0:
		raise ValueError('no pdbsolv file was found')

	return exposedhphob_path[0],bioplib_data_path[0],pdbsolv_path[0]

def get_happy(file):
	(exposedhphob_path,bioplib_data_path,pdbsolv_path)=get_file_paths()
	columns=["Happiness_mean","Nr_sad"]
	nr_sad=0
	#os.popen("export DATADIR={}".format(str(os.path.abspath("libs/bioplib-master/data"))))
	cmd="{} H95 H102 {}".format(exposedhphob_path,str(os.path.join(file)))
	command = subprocess.check_output(exposedhphob_path+" H95"+" H102 "+str(os.path.join(file)), shell=True, 
		env={'DATADIR': "{}".format(bioplib_data_path), "PDBSOLVDIR": pdbsolv_path},universal_newlines=True)
	command=command.splitlines()
	print(command)

	for i in command:
		print(i)
		i=i.split()
		print(i)
		if "Mean:" in i[0]:
			mean=i[1]
			#print(mean)
			continue
		if "Total:" in i[0]:
			pass
		elif float(i[2]) <0.5:
			nr_sad+=1
	try:
		write=[[mean,nr_sad]]
		happy_df = pd.DataFrame(write, columns=columns, index=None)
		print(happy_df)
		return (happy_df)
	except:
		print("error: happiness")
		return pd.DataFrame()
