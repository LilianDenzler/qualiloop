#!/usr/bin/env python
import sys
import os, shutil
import numpy as np
import csv
import pandas as pd
import matplotlib.pyplot as plt

from bs4 import BeautifulSoup
from urllib.request import urlopen

def write_alignment_input(sequence, name, input_seq_directory, file):
	full_seq=""
	with open(os.path.join(input_seq_directory,file)) as f:
		for line in f:
			line2=line.split()
			aminoacid=line2[1]
			full_seq+=str(aminoacid)
		with open(os.path.join("alignment_files","alignment_input_"+sequence+".txt"), 'a') as file:
			file.write(">{}\n".format(name))
			file.write("{}\n".format(full_seq))

def remove_small_files():
	for file in os.listdir(os.path.join("alignment_files")):
		file_path=os.path.join("alignment_files",file)
		file = open(file_path,"r")
		content = file.read()
		nr_lines = len(content.split("\n"))
		print("LENGTH:", nr_lines)
		if nr_lines<4: #as because of the newline also files with only only structure will technically have three lines
			print("REMOVED")
			os.remove(file_path)
		
def loop_redundancy(input_seq_directory,RMSD_file, save_path):
	if os.path.exists("alignment_files"):
		shutil.rmtree("alignment_files")
	os.makedirs(os.path.join("alignment_files"))
	list_seq={}
	doubles={}
	for file in os.listdir(input_seq_directory):
		with open(os.path.join(input_seq_directory,file)) as f:
			name=file.replace(".seq","")
			copy = False
			sequence=""
			for line in f:
				if "H95" in line:
					copy = True
				if copy==True:
					line2=line.split()
					aminoacid=line2[1]
					sequence+=str(aminoacid)
				if "H102" in line:
					copy = False
			if sequence in doubles.keys():
				doubles[sequence].extend([name])
				write_alignment_input(sequence, name,  input_seq_directory, file)

			elif sequence in list_seq:
				doubles.update({sequence:[name]})
				write_alignment_input(sequence, name,  input_seq_directory, file)

			else:
				list_seq.update({sequence:name})
				write_alignment_input(sequence, name,  input_seq_directory, file)

		f.close()
	RMSD_df=pd.read_csv(RMSD_file)
	range_list={}
	for key, values in doubles.items():
		if len(values)<=1:
			continue
		print("SEQ", key, values)
		RMSD_list=[]

		row=RMSD_df.loc[RMSD_df['seq'] == key]
		for index, i in row.iterrows():
			val=i["local_CA"]
			RMSD_list.append(val)
		print("RMSD_list", RMSD_list)
		if len(RMSD_list)<2:
			continue
		range_RMSD=max(RMSD_list)-min(RMSD_list)
		print("RANGE",range_RMSD)
		range_list[range_RMSD]=[key, values]
	max_key=max(range_list.keys())
	print("MAAAX",range_list[max_key])
	print("MAAAX",max_key)
	y=list(range_list.keys())
	X=[i[0] for i in range_list.values()]

	plt.plot(X, y, 'o', color='black')
	plt.ylabel('Range of RMSD values')
	plt.xlabel("Loop Sequence")
	plt.xticks(rotation=90, fontsize=10)
	plt.title('RMSD Ranges of Duplicated Loop Sequences')
	plt.tight_layout()
	plt.savefig(save_path)
	remove_small_files()

def check_empty(to_be_checked_file):
	file = pd.read_csv(to_be_checked_file)
	with open(to_be_checked_file,'r') as csvfile: 
		reader = csv.reader(csvfile)
		counter=0

		for row in reader:
			if all(row[0:file.shape[1]]) and len(row)==file.shape[1]:
				pass
			else:
				counter+=1
	csvfile.close()

	print("From ",file.shape[0],"rows, " ,counter," rows will be deleted")
	file=file.dropna()
	file.to_csv(to_be_checked_file+".nonan",index=False)




def check_redundancy(to_be_checked_file, input_seq_directory):
	list_doubles1, list_doubles2= loop_redundancy(input_seq_directory)
	lines=[]
	print(to_be_checked_file)
	data = pd.read_csv(to_be_checked_file, error_bad_lines=False, header=0)
	x=data.target.tolist()
	for i in x:
		lines.append(i)
	for line in lines:
		for i in list_doubles1:
			if str(i) in line:                                            
				print (line)
				print(str(i))
				lines.remove(line)
				#print("hey")

	with open(to_be_checked_file+".nonan",'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(lines)
		writeFile.close()

def check_too_short(to_be_checked_file):
	data = pd.read_csv(to_be_checked_file, error_bad_lines=False, header=0)
	print("before shape: ", data.shape)
	print(data.length.tolist())
	data = data[data.length  != 1]
	data = data[data.length  != 2]
	data.to_csv(os.path.join(to_be_checked_file+".noshort"),index=False)
	print("after shape: ", data.shape)
   

def duplicate_IDs(to_be_checked_file):
	data = pd.read_csv(to_be_checked_file)
	data.drop_duplicates(subset ="ID", keep = False, inplace = True) 
	print(data.shape)
	print(data)
	data.to_csv(to_be_checked_file+".fixed",index=False)



def resolution_scraper(queries):
	for query in queries:
		url = 'https://www.rcsb.org/structure/' + querys
		html=urlopen(url)
		bsObj=BeautifulSoup(html)
		
		resolution=bsObj.findAll(id="exp_header_0_diffraction_resolution")
		resolution=resolution[0].get_text().replace("&nbsp", "")

		r_value=bsObj.findAll(id="exp_header_0_diffraction_rvalueFree")
		r_value=r_value[0].get_text().replace("&nbsp", "")
		print(resolution)
		print(r_value)
		entity_dic={}
		for i in range(5):
			try:
				table=bsObj.find("table",{"id":"table_macromolecule-protein-entityId-{}".format(i)})
				for row in table.find_all('tr', {"id": "macromolecule-entityId-{}-rowDescription".format(i)}):
					columns = row.find_all('td')
					entity_list=[]
					for a in columns:
						add=a.get_text().split("&nbsp")[0]
						entity_list.extend([add])

					entity_list=entity_list[:-1]
					entity_dic[i]=entity_list
			except:
				pass
		entity_df=pd.DataFrame.from_dict(entity_dic, columns=["Molecule","Chains","Length","Organism",'Mutations'],orient='index')
		print(entity_df)

if __name__ == '__main__':
	#check_redundancy(sys.argv[1], sys.argv[2])
	#check_empty(sys.argv[1])
	#duplicate_IDs(sys.argv[1])
	loop_redundancy(sys.argv[1], sys.argv[2], sys.argv[3])
	#check_too_short(sys.argv[1])
	#resolution_scraper()

