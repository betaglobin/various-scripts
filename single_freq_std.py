#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 13:57:11 2018

@author: staskowiak
"""
import math


AA_dict = { "A" : "Ala", "R" : "Arg", "D" : "Asp", "N" : "Asn", "C" : "Cys",
            "E" : "Glu", "Q" : "Gln", "G" : "Gly", "H" : "His", "I" : "Ile", 
            "L" : "Leu", "K" : "Lys", "M" : "Met", "F" : "Phe", "P" : "Pro", 
            "S" : "Ser", "T" : "Thr", "W" : "Trp", "Y" : "Tyr", "V" : "Val"}


glob_freq = {"A": 8.25, "R": 5.53, "D": 5.45, "N": 4.06, "C": 1.37,  
             "E": 6.75, "Q": 3.93, "G": 7.07, "H": 2.27, "I": 5.96,
             "L": 9.66, "K": 5.84, "M": 2.42, "F": 3.86, "P": 4.70, 
             "S": 6.56, "T": 5.34, "W": 1.08, "Y": 2.92, "V": 6.87}



file_name = "Arabidopsis_thaliana.TAIR10.pep.all.fa"
fp = open(file_name, "r")
full_file = fp.read()

full_heads = []
full_seqs = []

ctrl = 0
i = 0
while True:
    head = ""
    while True:
        try:
            head = head + full_file[i]
            i +=1
            if full_file[i] == "\n":
                break
            
        except IndexError:
            break
            
    seq  = ""
    while True:
        try:
            seq = seq + full_file[i]
            i += 1
            if full_file[i] == ">":
                break
    
        except IndexError:
            ctrl = 1
            break
        
    if ctrl == 1:
        break

    full_heads.append(head)
    full_seqs.append(seq)
    
    
for i in range(len(full_seqs)):
    seq = ""
    for j in full_seqs[i]:
        if j != "\n":
            seq = seq + j
            
    full_seqs[i] = seq



full_seqs


AA_array = ["A", "R", "D", "N", "C", "E", "Q", "G",
            "H", "I", "L", "K", "M", "F", "P", "S",
            "T", "W", "Y", "V"]


def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

numerator = []
for i in range(len(AA_array)):
    numerator.append(0)
    
    
denominator = len(find(full_file, ">"))



for i in range(len(AA_array)):  
    
    AA_occur = []
    for j in range(len(AA_array)):
        AA_occur.append(0)          #prepering array
    
    for ch in full_seqs[i]:
        AA_occur[AA_array.index(ch)] += 1   #update array
        

    l = 0
    for k in AA_occur:
        substract = k-glob_freq[AA_array[l]] 
        result = math.pow(substract, 2)    #add it to the numerator of the variance
        numerator[l] += result
        l += 1
        
   
st_dev = []

for i in numerator:
    st_dev.append(round(math.sqrt(i/denominator),2))
    
    


#saving stuff to cvs
row = "amino_acid;mean_freq;std\n"
fp = open("AA_content.txt", "a")
fp.write(row)  
  
   
for i in range(len(AA_array)):
    row = ""
    fp = open("AA_content.txt", "a")
    row = AA_array[i] + ";" + str(glob_freq[AA_array[i]]) + ";" + str(st_dev[i]) + "\n"
    fp.write(row)








