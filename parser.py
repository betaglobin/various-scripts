#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 14:10:16 2018

@author: staskowiak
"""
import os

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

#extract two gene names from the headers and check if those are the same
def compare(x, y):
    header1 = ""
    j = x    
    while(True):
        header1 = header1 + full_file[j]
        if (full_file[j] == "\n"):
            break
        j += 1
              
    k = header1.find("gene:") + 5
    gname1 = ""
    
    while(True):
        gname1 = gname1 + header1[k]
        k += 1
        if (header1[k] == " "):
            break
    
    if y != "abc":
        header2 = ""
        j = y
        while(True):
            header2 = header2 + full_file[j]
            if (full_file[j] == "\n"):
                break
            j += 1
                    
        k = header2.find("gene:") + 5
        gname2 = ""
        
        while(True):
            gname2 = gname2 + header2[k]
            k += 1
            if (header2[k] == " "):
                break
                 
        #compare gene names       
        if gname1 == gname2:
            return True
        
        if gname1 != gname2:
            return False
        
    elif y == "abc":
        return False
    
    
global full_file
full_file = ""
fp = open("roboczy.fa", "r")
full_file = fp.read()
fp.close()

positions = find(full_file, ">") #establish the positions of the sequences

parsed_array = []
queue_temp = []
i = 0
extend = 1
while True:  
    try:
        fst = positions[i]        #get the indexes...
        try:
            scd = positions[i+extend]
        except IndexError:
            scd = "abc" #a control variable if you reached EOF
        
        result = compare(fst, scd) #... and make a comparison
        
        queue_temp.append(fst)
        
        if result == True:
            queue_temp.append(scd) #if the names are equal, create a family of isoforms
            i += 1

        if result == False: 
            queue_temp = list(set(queue_temp)) 
            parsed_array.append(queue_temp) #if not then add this family to the array
            i += 1
            extend = 1     
            queue_temp = []
            
    except IndexError:
        break
        
    
def seq_extractor(x):
    while(True):                    #starting from the header...
        if full_file[x] != "\n":
            x += 1
        else:
            break
        
    seq_beg = x + 1                 #... establish postion of the sequence
    f_seq = ""
    
    while(True):
        if full_file[seq_beg] != ">" or full_file[seq_beg] != "\n":
            f_seq = f_seq + full_file[seq_beg]  #don't take blank characters
            
        seq_beg += 1 
        
        try:
            if full_file[seq_beg] == ">":
                break
        except IndexError:
            break
        
    return f_seq
            
                  
#parsed_array

uniq_array = []

for i in parsed_array:
    if len(i) == 1:
        uniq_array.append(i)
    else:
        keys = []
        vals = []
        for j in i:
            keys.append(j)
            seq = seq_extractor(j)
            vals.append(seq)
            
        x = 0
        ind = 0
        for k in range(len(vals)):
            curr = len(vals[k])
            if (curr > x):
                x = curr
                ind = k
    
        temp = []
        temp.append(keys[ind])
        uniq_array.append(temp)
                                       
#uniq_array
#len(uniq_array)

for i in uniq_array:
    gene = ""
    k = i[0]
    
    while(True):
        gene = gene + full_file[k]
        k += 1
        
        try:
            if (full_file[k] == ">"):
                break
        except IndexError:
            break
    
    fp = open("roboczy.FILTERED.fa", "a")
    fp.write(gene)
    fp.close()

