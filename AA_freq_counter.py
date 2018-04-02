#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 18 17:40:04 2018

@author: staskowiak
"""

AA_array = ["A", "R", "D", "N", "C", "E", "Q", "G",
            "H", "I", "L", "K", "M", "F", "P", "S",
            "T", "W", "Y", "V"]

global full_file
full_file = ""
file_name = input("Provide a filename:")


fp = open(file_name, "r")
full_file = fp.read()
fp.close()


#read everything to two lists: one containing headers and the other containing sequences
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
    

#removing blank characters from the sequences
for i in range(len(full_seqs)):
    seq = ""
    for j in full_seqs[i]:
        if j != "\n":
            seq = seq + j
            
    full_seqs[i] = seq


#create a list containing gene ID in Uniprot format
gene_array = []

for head in full_heads:
    gstart = head.find("gene:") + 5
    
    gname = ""
    while True:
        if head[gstart] == " ":
            break
        
        gname = gname + head[gstart]
        gstart += 1

    gene_array.append(gname)


#extracting gene symbol (gene name) if exists
gene_symbol = []
for i in full_heads:
    pos = i.find("gene_symbol")
    
    g_name = ""
    if pos == -1:
        g_name = "NULL"
        
    if pos > -1:
        pos += 12
        k = pos
        while True:
            try:
                if i[k] == " ":
                    break
                k += 1
                
            except IndexError:
                break
            
        g_name = i[pos:k]
        
    gene_symbol.append(g_name)
        



AA_pairs = []

i = 0
j = 0

while True:    
    while True:
        pair = ""
        try:
            pair = AA_array[j] + AA_array[i]
            i += 1
            AA_pairs.append(pair)
        except IndexError:
            break
        
    i = 0
    j += 1
        
    if j == len(AA_pairs):
        break
     



#making a header of a CVS file
AA_feature = ""
for i in AA_pairs:
    AA_feature = AA_feature + i + ";"

AA_feature = AA_feature[:-1]
header = "gene;gene_symbol;sequence;" + AA_feature
header += "\n"

fp = open("data_freq_pairs_COMPLETE.txt", "w")
fp.write(header)


l = 0


for seq in full_seqs:
    single_rec = ""
    single_rec = gene_array[l] + ";" + gene_symbol[l] + ";" + seq + ";"
    
    for pair in AA_pairs:
        k = 0
        AA_count = 0
        
        while True:
            try:
                cmpr = seq[k]+seq[k+1]
                if cmpr == pair:
                    AA_count +=1

                k +=1
            except IndexError:
                break
            
        single_rec = single_rec + str((AA_count/(len(seq)-1))) + ";" 
    
    single_rec = single_rec[:-1]    
    single_rec += "\n"
    
    fp = open("data_freq_pairs_COMPLETE.txt", "a")
    fp.write(single_rec)
    l +=1

