#!/usr/bin/env python
"""/*************************************************************************
   Program:    Sequence Parser
   File:       seq_parser.py
   
   Version:    V0.2
   Date:       10.10.20
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
   The program takes a directory of experimentally determined PDB structures 
   and outputs sequence files foe each of the input PDB files. The generated 
   sequence files are to be used as input files for abYmod. 
**************************************************************************
   Usage:
   ======
   This program is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   10.10.20 Original
   V0.2	  09.05.21 Revision
*************************************************************************/"""

import os
import sys

light=[]
heavy=[]

#argv[1]= direcoty with actual PDBs, argv[2] = directory to keep input sequences

def one_letter_code(residue):
	dic = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M','XAA': 'X', 'UNK':'X'}
	#one_letter= atomium.structures.code()    #####doesnt work -> why?
	if len(residue) % 3 != 0:
		print(residue)
		raise ValueError("error")
	one_letter= dic[residue]
	return (one_letter)

def parser(position,residue,chain_type):
	global string
	for i in zip(position, residue, chain_type):
		one_res=one_letter_code(residue)
		if chain_type=="H" or chain_type=="L":
			string+= '{}{}\t{}\n'.format(chain_type,position,one_res)
	return (string)

if __name__ == '__main__':
	for filename in os.listdir(sys.argv[1]):
		string=[]
		file=open(sys.argv[1]+"/"+filename)
		if filename.endswith(".pdb")== False:
			print("this file is not in correct format (i.e. not PDB){}".format(filename))
			continue
		for line in file:
			if "ATOM" in line and "CA" in line:
				fields = line.strip().split()
				pos= fields[5]
				res= fields[3]
				chain= fields [4]
				string= parser(pos, res, chain)
		file.close()

		new_input_file= open("{}/{}.seq".format(sys.argv[2], filename[:-4]),"w")
		for i in string:
			new_input_file.write(i)
