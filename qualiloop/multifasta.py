#!/usr/bin/python3
"""/*************************************************************************
   Program:    Multi-Fasta Maker
   File:       multifasta.py
   
   Version:    V1.1
   Date:       30.06.21
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
   The functions contained here are to be used to create multiple-Fasta
   file from directories containing a number of PDB structures.
   The program was written to enable abYmod vs PIGS comparison of 
   CDRH3-structure modelling.
**************************************************************************
   Usage:
   ======
   This library is intended for the Qualiloop application. 
**************************************************************************
   Revision History:
   =================
   V0.1   30.06.21 Original
*************************************************************************/
example output:

>Immunoglobulin_1 VL
DIVLTQSPASLSASVGETVTITCRASGNIHNYLAWYQQKQGKSPQLLVYYTTTLADGVPS
RFSGSGSGTQYSLKINSLQPEDFGSYYCQHFWSTPRTFGGGTKLEIK
>Immunoglobulin_1 VH
QVQLQESGPGLVAPSQSLSITCTVSGFSLTGYGVNWVRQPPGKGLEWLGMIWGDGNTDYN
SALKSRLSISKDNSKSQVFLKMNSLHTDDTARYYCARERDYRLDYWGQGTTLTVSS
"""
import sys
import os
import pandas as pd
import numpy as np

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

def get_heavy_light(actual_dir, run_dir):
   open(os.path.join(run_dir,"multifasta.txt"), 'w')
   for file in os.listdir(actual_dir):
      file_name=file.split(".")[0]
      file=os.path.join(actual_dir, file)
      input_seq_df=get_df(file)
      light_chain=input_seq_df[input_seq_df['chain'].astype(str) == "L"]["res"].to_list()
      heavy_chain=input_seq_df[input_seq_df['chain'].astype(str) == "H"]["res"].to_list()
      with open(os.path.join(run_dir,"multifasta.txt"), 'a') as file:
         if len(heavy_chain)==0 or len(light_chain)==0:
            pass
         else:
            file.write(">"+file_name+" VL \n")
            file.write("".join(light_chain)+"\n")
            file.write(">"+file_name+" VH \n")
            file.write("".join(heavy_chain)+"\n")
   return (input_seq_df)


if __name__ == '__main__':
   get_heavy_light(sys.argv[1], sys.argv[2])