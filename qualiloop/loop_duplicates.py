#!/usr/bin/env python
import sys
import os
import numpy as np
import csv
import pandas as pd

"""loop through all lines of the sequence file, extract the loop sequence (H95-H102)
            each sequence file has the following format:
            L1  N
            L3  V
            L4  L
            ...
            H100  N
            ...
            the goal is to extract the single-letter AA code for residues H95-H102.
            Then, these sequences are stores in a dictionary. If a sequence is already found in the dictionary, it will be added to the doubles dictionary.
            This is done to identify the redundant loop sequences in the dataset.
            The doubled sequences are then 
            """   
def loop_redundancy(input_seq_directory):
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
            if sequence in list_seq:
                doubles.update({sequence:name})
            else:
                list_seq.update({sequence:name})
        f.close()
    for keys,values in doubles.items():
        doubles[keys]=[doubles[keys],list_seq[keys]]
    list_doubles1=[]
    list_doubles2=[]
    for keys, values in doubles.items():
        b=doubles[keys]
        a=b[0]
        list_doubles1.append(a)
    for keys, values in doubles.items():
        b=doubles[keys]
        a=b[1]
        list_doubles2.append(a)
    for keys, values in doubles.items():
        b=doubles[keys]
        additionals=[]
        if len(b)>2:
            additionals+=b[2,:]
    return (list_doubles1, list_doubles2, additionals)

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
    list_doubles1, list_doubles2, additionals= loop_redundancy(input_seq_directory)
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






if __name__ == '__main__':
    check_empty(sys.argv[1])
