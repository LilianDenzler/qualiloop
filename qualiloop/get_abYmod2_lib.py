#!/usr/bin/python3
#python3 get_abYmod2.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/
#nohup python3 get_abYmod2.py /home/lilian/sync_project/Redundant_PDBs.txt /home/lilian/sync_project/input_Abymod_seq/ /home/lilian/sync_project/abYmod_structures/ &>output_abymod.log &
import sys
import os
import shlex
#import subprocess
from subprocess import Popen, PIPE, call
import csv


def pass_commands(seq_csv, upload_dir, ID, Redundant_PDBs=None):
	with open(Redundant_PDBs, "r") as f:
	    reader = csv.reader(f)
	    redundants = {}
	    for row in reader:
	        key = row[0]
	        redundants[key] = row[1:]
	to_exclude=str(redundants.get(ID)).strip('[]').replace("'","").replace(" ", "").rstrip(',')
	
	#os.system("abymod -v=4 -excl100 {} > {} 2>{}".format(seq_csv, os.path.join(upload_dir,ID+".pdb.model"), os.path.join(upload_dir,ID+".log")))
	#print("abymod -v=4 -excl100 {} > {} 2>{}".format(seq_csv, os.path.join(upload_dir,ID+".pdb.model"), os.path.join(upload_dir,ID+".log")))
	command="abymod -v=4 -autoexclude {} > {} 2>{}".format(to_exclude, seq_csv, str(os.path.join(upload_dir,ID+".pdb.model")), str(os.path.join(upload_dir,ID+".log")))
	#-excl100 -exclude={}
	#command=shlex.split(command)
	#csv=os.path.join(upload_dir,ID+".abymod_out")
	print(command)
	p=Popen(command, shell=True)
	p.wait()
	return (os.path.join(upload_dir,ID+".pdb.model"),os.path.join(upload_dir,ID+".log")) 
