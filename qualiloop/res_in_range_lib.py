import sys
sys.path.append('/serv/www/html_lilian/libs')
sys.path.append('./CDRH3lib')
import os
import pandas as pd
import numpy as np
import subprocess
import glob 

def get_file_paths():
	rangecontacts_path=glob.glob('**/libs/**/rangecontacts', recursive=True)
	if len(rangecontacts_path)>1:
		print(rangecontacts_path)
		#raise Warning('more than one range contacts file was found')
	elif len(rangecontacts_path)==0:
		print("rangecontacts file is not in lib!!!!")
		rangecontacts_path=glob.glob('**/rangecontacts', recursive=True)
		raise Warning("rangecontacts file is not in lib!!!!")
		if len(rangecontacts_path)==0:
			raise ValueError('no rangecontacts file was found')
	return rangecontacts_path[0]


def get_in_range(file):
	rangecontacts_path=get_file_paths()
	#get the number of contacts the loop makes within and outside of the loop
	print(rangecontacts_path)
	command=[rangecontacts_path, "-i", "-m", "H95","H102",file]
	proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	outputs_raw, errs = proc.communicate()
	print(outputs_raw, errs)
	rangecontacts_in=outputs_raw.splitlines()

	command=[rangecontacts_path, "-m", "H95","H102",file]
	proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	outputs_raw, errs = proc.communicate()
	print(outputs_raw, errs)
	rangecontacts_out=outputs_raw.splitlines()


	#command = subprocess.check_output(str(rangecontacts_path+ " -i"+ " -m"+ " H95"+" H102"+file), shell=True, universal_newlines=True)
	#rangecontacts_in=command.splitlines()
	#print(command)

	#get the number of contacts the loop makes only outside the loop
	#command = subprocess.check_output(rangecontacts_path+" -m"+ " H95"+" H102"+file, shell=True, universal_newlines=True)
	rangecontacts_in=outputs_raw.splitlines()
	in_range_df=pd.DataFrame()
	if len(rangecontacts_in)==0:
		return pd.DataFrame()
	in_range_df["contacts_all"]=[len(rangecontacts_in)]
	in_range_df["contacts_out"]=[len(rangecontacts_out)]
	in_range_df["contacts_ratio_out_all"]=[len(rangecontacts_out)/len(rangecontacts_in)]
	return in_range_df


if __name__ == '__main__':
	in_range_df=get_in_range(sys.argv[1])
	print(in_range_df)