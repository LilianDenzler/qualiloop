#!/usr/bin/python3
import sys
import os
import numpy as np
import csv
import pandas as pd
import tempfile
import subprocess

def get_avp_input(file):
	file=open(file)
	copy=False
	seq=[]
	for line in file:

		if "ATOM" in line:
			if "H  95" in line:
				copy = True
				seq.append(line)
				continue
			elif "H 103" in line:
				#seq.append(line)
				copy = False
				continue
			elif copy==True:
				seq.append(line)
	loop_pdb_tmp = open("./loop_pdb_tmp.txt", "w")
	for i in seq:
		loop_pdb_tmp.write(i)
	#print(loop_pdb_tmp.read())
	loop_pdb_tmp.close()
	return None
		
	
def run_avp():
	#make temporary files to store the avp output
	avp_voids_tmp = open("./avp_voids_tmp.txt", "a")
	near_voids_tmp = open("./near_voids_tmp.txt", "a")
	avp_voids_tmp.close()
	near_voids_tmp.close()

	#run the avp command - -n option so that additional file is outputted that contains the pdb data of all atoms near the voids
	#solvent size of 1.4 Angstroms (water-molecule) and probe size of 0.5 is suggested in the paper by A. CR Martin (as if probe too small, 
	#separate voids join into one huge one)
	command="avp -q -R -p 0.5 -s 1.4 -n {} {} {}".format("./near_voids_tmp.txt", "./loop_pdb_tmp.txt", "./avp_voids_tmp.txt")
	avp_out_tmp=subprocess.call(command, shell=True)

	#read and parse the avp output file with the voids
	columns_voids=["largest_total", "buried_surface_void", "type", "overall_vols", "surface_buried_vols", "ind_vols", "midpoint", "x", "y", "z", "void_type"]
	voids=pd.read_csv("./avp_voids_tmp.txt", header=None, delim_whitespace=True, names=columns_voids)
	if len(voids.loc[(voids['largest_total'] == "Void:")].to_numpy())==0:
		return pd.DataFrame()
	surface_ind_voids= voids.loc[(voids['midpoint'] == "midpoint")&(voids['void_type']=="Buried-void")]
	surface_ind_voids=surface_ind_voids[["x", "y", "z", "ind_vols"]]
	buried_ind_voids= voids.loc[(voids['midpoint'] == "midpoint")&(voids['void_type']=="Buried-void")]
	buried_ind_voids=buried_ind_voids[["x", "y", "z", "ind_vols"]]
	summarized_voids=pd.DataFrame()

	largest_void=(voids.loc[(voids['largest_total'] == "Largest")&(voids['buried_surface_void']=="void"),"overall_vols"]).to_numpy()[0]
	summarized_voids["total_buried"]=voids.loc[(voids['largest_total'] == "Total")&(voids['buried_surface_void']=="buried"), "surface_buried_vols"]
	summarized_voids["total_surface"]=float(voids.loc[(voids['largest_total'] == "Total")&(voids['buried_surface_void']=="surface"),"surface_buried_vols"])
	summarized_voids["largest_buried"]=float(voids.loc[(voids['largest_total'] == "Largest")&(voids['buried_surface_void']=="buried"),"surface_buried_vols"])
	summarized_voids["largest_surface"]=float(voids.loc[(voids['largest_total'] == "Largest")&(voids['buried_surface_void']=="surface"),"surface_buried_vols"])
	summarized_voids["total_void"]=float(voids.loc[(voids['largest_total'] == "Total")&(voids['buried_surface_void']=="void"),"overall_vols"])
	summarized_voids["largest_void"]=float(voids.loc[(voids['largest_total'] == "Largest")&(voids['buried_surface_void']=="void"),"overall_vols"])
	largestx=voids.loc[(voids['largest_total'] == "Void:")&(voids['ind_vols']==max(voids["ind_vols"].to_numpy()))]
	summarized_voids["largest_x"]=largestx[["x"]].to_numpy()
	largesty=voids.loc[(voids['largest_total'] == "Void:")&(voids['ind_vols']==max(voids["ind_vols"].to_numpy()))]
	summarized_voids["largest_y"]=largesty[["y"]].to_numpy()
	largestz=voids.loc[(voids['largest_total'] == "Void:")&(voids['ind_vols']==max(voids["ind_vols"].to_numpy()))]
	summarized_voids["largest_z"]=largestx[["z"]].to_numpy()

	#read and parse the avp output file containing the atoms near the voids
	columns_near=["atom_nr", "atom_name", "res", "chain", "pos","x", "y", "z", "occ", "temp", "element"]
	near_voids=pd.read_csv("./near_voids_tmp.txt", header=None, delim_whitespace=True, names=columns_near)
	summarized_near=pd.DataFrame()
	summarized_near=near_voids[["res", "atom_name","chain", "pos", "x", "y", "z"]]
	dist_list=[]
	for i in range(summarized_near.shape[0]):
		nearest=(abs(float(summarized_near["x"].to_numpy()[i])-float(summarized_voids["largest_x"].to_numpy()[0]))
		+abs(float(summarized_near["y"].to_numpy()[i])-float(summarized_voids["largest_y"].to_numpy()[0]))
		+abs(float(summarized_near["z"].to_numpy()[i])-float(summarized_voids["largest_z"].to_numpy()[0])))
		
		dist_list.append(nearest)
	index_min=dist_list.index(min(dist_list))
	closest=summarized_near.iloc[index_min].to_numpy()
	closest = pd.DataFrame(closest.reshape(-1, len(closest)),columns=["res", "atom_name","chain", "pos", "x_closest", "y_closest", "z_closest"], index=None)
	print(closest)
	
	#summarized_near=summarized_near.drop_duplicates(subset=["res", "chain", "pos"])

	#join the datafram on the voids with the df on the atom closest to the largest void

	summarized_voids.reset_index(drop=True, inplace=True)
	closest.reset_index(drop=True, inplace=True)
	voids_df=pd.concat([summarized_voids, closest], axis=1)
	

	#remove the files, as no longer needed
	os.remove("./near_voids_tmp.txt")
	os.remove("./avp_voids_tmp.txt")
	os.remove("./loop_pdb_tmp.txt")

	return voids_df

def get_voids(file):
	get_avp_input(file)
	voids_df=run_avp()
	voids_df=voids_df[["pos", 'largest_void', "total_surface", "largest_surface", "total_void"]]
	print(voids_df)
	return voids_df


if __name__ == '__main__':
	print(get_voids(sys.argv[1]))



#sample output from avp for the voids
"""Void:    1 voxelCount:   3 volume:    2.218 midpoint:   21.825  40.522  31.265 Surface-void
Void:    2 voxelCount:   1 volume:    2.342 midpoint:   23.825  40.022  32.765 Surface-void
Void:    3 voxelCount:  12 volume:    9.664 midpoint:   26.325  40.022  30.765 Surface-void
Void:    4 voxelCount:   1 volume:    0.922 midpoint:   24.825  41.022  32.765 Surface-void

Total buried void volume: 0.000000
Largest buried void volume: 0.000000

Total surface void volume: 15.146000
Largest surface void volume: 9.664000

Total void volume: 15.146000
Largest void volume: 9.664000"""

#sample output from avp for the nearest atoms to voids
"""ATOM   1589  O   ARG H  96      23.300  40.621  30.807  1.00  0.00           O  
ATOM   1597  O   ASP H  97      21.354  39.948  33.284  1.00  0.00           O  
ATOM   1610  N   ARG H  99      24.506  38.512  33.129  1.00  0.00           N  
ATOM   1620  O   ARG H  99      26.261  38.955  30.956  1.00  0.00           O  
ATOM   1589  O   ARG H  96      23.300  40.621  30.807  1.00  0.00           O  
ATOM   1632  CG  ASP H 101      26.201  39.441  26.525  1.00  0.00           C  
ATOM   1571  CA  GLU H  95      28.096  41.866  31.323  1.00  0.00           C  
ATOM   1610  N   ARG H  99      24.506  38.512  33.129  1.00  0.00           N  
ATOM   1579  N   ARG H  96      25.932  41.924  30.352  1.00  0.00           N  
ATOM   1598  N   TYR H  98      23.252  40.828  34.152  1.00  0.00           N"""