#!/usr/bin/python3
import sys
import os
import pandas as pd
import subprocess
import glob


def get_file_paths():
	pdbsolv_path=glob.glob('**/libs/**/bioptools-master/src/pdbsolv', recursive=True)
		
	if len(pdbsolv_path)>1:
		print(pdbsolv_path)
		#raise Warning('more than one pdbsolv file was found')
	elif len(pdbsolv_path)==0:
		print("pdbsolv file is not in lib!!!!")
		rangecontacts_path=glob.glob('**/pdbsolv', recursive=True)
		#raise Warning("pdbsolv file is not in lib!!!!")
		if len(pdbsolv_path)==0:
			raise ValueError('no pdbsolv file was found')

	
	bioplib_data_path=glob.glob('**/libs/**/bioplib-master/data', recursive=True)
	if len(bioplib_data_path)>1:
		raise ValueError('more than one bioplib data file was found')
	elif len(bioplib_data_path)==0:
		print("bioplib data is not in lib!!!!")
		rangecontacts_path=glob.glob('**/bioplib-master/data', recursive=True)
		#raise Warning("bioplib data is not in lib!!!!")
		if len(bioplib_data_path)==0:
			raise ValueError('no bioplib data file was found')

	return pdbsolv_path[0], bioplib_data_path[0]

def pdbsolv_run(file):
	columns=["Access","Relacc","Scacc","Screlacc","Access_avg","Relacc_avg","Scacc_avg","Screlacc_avg"]
	data=[]
	filename=file
	nr_sad=0
	pdbsolv_path, bioplib_data_path=get_file_paths()

	print(pdbsolv_path+" -r"+" stdout "+"{}".format(os.path.join(file)))
	command = subprocess.check_output(pdbsolv_path+" -r"+" stdout "+"{}".format(os.path.join(file)), shell=True, 
		env={'DATADIR': "{}".format(bioplib_data_path), "PDBSOLVDIR": pdbsolv_path},universal_newlines=True)

	command=command.splitlines()
	#cmd_in="pdbsolv -r stdout {}".format(os.path.join(file))
	#command = subprocess.Popen(cmd_in, shell=True, env={'DATADIR': '/serv/www/html_lilian/libs/data'}).readlines()
	# RESIDUE  AA   ACCESS  RELACC  SCACC   SCRELACC
	copy=False
	results=False
	access=0
	relacc=0
	scacc=0
	screlacc=0
	counter=0
	for i in command:
		i=i.split()
		#print(i)
		if "END" in i:
			results=True
			continue
		if results==True:
			if "95"== i[2] and "H"==i[1]:
				copy = True
			if copy==True:
				counter+=1
				access+=float(i[4])
				relacc+=float(i[5])
				scacc+=float(i[6])
				screlacc+=float(i[7])
			if "102"==i[2] and "H"==i[1]:
				copy = False
	if counter==0:
		access_avg=None
		relacc_avg=None
		scacc_avg=None
		screlacc_avg=None
		access=None
		relacc=None
		scacc=None
		screlacc=None
		write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg]
		data.append(write)

	else:
		access_avg=access/counter
		relacc_avg=relacc/counter
		scacc_avg=scacc/counter
		screlacc_avg=screlacc/counter

		write=[access,relacc,scacc,screlacc, access_avg,relacc_avg, scacc_avg, screlacc_avg]
		data.append(write)
	
			
	df = pd.DataFrame(data, columns=columns)
	return (df)

if __name__ == '__main__':
	df=pdbsolv_run(os.path.join("/serv/www/html_lilian","1F11_1.pdb"))
	print(df)