#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: staskowiak

#the following script parses a file at:
#https://raw.githubusercontent.com/aziele/bsb-2018/master/data/kinases.genbank.txt
#gets headers and sequences and saves those in FASTA format

full_file = ""
full_file = open("data.txt", "r")
full_file = full_file.read()

backup = full_file

pre_array = []
pre_concat = ""


#       1. SEQUENCE EXTRACTION

#find a position of the ORIGIN in order to get to the sequence
#after that trimm the file and reapet
while(full_file.find("ORIGIN") > 0):   
    index = full_file.find("ORIGIN")
    if index == -1:
        break
    index = index + 8 #move fixed number of characters to the right
    pre_concat = ""
    while(full_file[index] != "/"):
        pre_concat = pre_concat + full_file[index]
        index = index+1
        
    pre_array.append(pre_concat)
    pre_concat = ""
    
    full_file = full_file[index:] #cut the file
    
pre_array

import string
seq_array = []

for i in pre_array:
    empty_seq = ""
    n = 0;
    while(True):
        if n == len(i):
            seq_array.append(empty_seq)
            empty_seq = ""
            break
    
        #i'm intersted in a-z characters only
        if i[n] in string.ascii_lowercase:
            empty_seq = empty_seq + i[n]
        
        n = n+1
             
#a conversion to "nice-to-have-in-fasta-file" capital letters              
for i in range(len(seq_array)):
    seq_array[i] = seq_array[i].upper()      

del(pre_array)
del(full_file)        
        

#                       2. HEADER EXTRACTION
headers = []
pre_concat = ""

while(backup.find("LOCUS") >= 0):
    index = backup.find("LOCUS")
    if index == -1:
        break
    
    pre_concat = ">"
    index = index + 12 #move fixed number of characters
    
    while(backup[index] != " "):
        pre_concat = pre_concat + backup[index]
        index = index+1
     
    pre_concat += "\n"
    headers.append(pre_concat)
    pre_concat = ""
    backup = backup[index:]
    
del(backup)
    
#headers

final_file = ""
index = 0
for i in headers:
    final_file +=i + seq_array[index] + "\n" + "\n"
    index += 1
    
    
f = open("parsed_results.fa", "a")
f.write(final_file)
f.close()
