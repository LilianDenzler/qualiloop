#!/usr/bin/env python3

"""
Program:    seq_parser_lib
File:       seq_parser_lib.py

Version:    V1.0
Date:       20.10.2020
Function:   Parse PDB sequence
            
Description:
============
This program takes a PDB file and extracts the CA atoms. It outputs a python pandas dataframe containing the chain, position and residue.
e.g.:
H 	1	L
H 	2	Q
...
--------------------------------------------------------------------------
Revision History:
=================
V1.0   20.10.2020
"""
#*************************************************************************
# Import libraries
import pandas as pd
import sys 
#*************************************************************************


def one_letter_code(residue):

	"""
	Go from the three-letter code to the one-letter code.

	Input:  residue      --- The three-code residue identifier
    Return: one_letter --- The one-letter residue identifier

    """

	dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA': 'X', 'UNK':'X'}
	if len(residue) % 3 != 0:
		if residue[1:] in dic.keys() and residue[0] in ["A","B", "C", "D"]:
			one_letter= dic[residue[1:]]
			return one_letter
		print(residue)
		raise ValueError("error")
	one_letter= dic[residue]
	return one_letter


def get_df(file):
	"""
	Take a PDB file and make a dataframe containing all CA-atoms with information on chain, position and residue

	Input:  PDB file     
    Return: dataframe  --- columns: chain, pos, res

    """
	dic_input= {'chain': [],'pos':[],'res':[]}
	with open(file, "r") as file:
		for line in file:
			if "ATOM" in line and "CA" in line:
				if "REMARK" in line:
					continue
				fields = line.strip().split()
				pos= fields[5]
				res= fields[3]
				chain= fields [4]
				one_res=one_letter_code(res)
				dic_input['chain'].append(chain)
				dic_input['pos'].append(pos)
				dic_input['res'].append(one_res)
	df = pd.DataFrame(dic_input, columns = ['chain', 'pos','res'], index=None)
	return (df)

if __name__ == '__main__':
	df=get_df(sys.argv[1])
	print(df)