#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 11:14:21 2023

Script that will read through construct files with values and count how often
each individual gene spacer (gene ID) is found in any combination of spacer
constructs. Written for Nicole Ellis for her eLife paper by Kevin Myers.

Requires argparse, re, sys

Input takes two files:
    1) spacer_ids.txt = a list of all gene IDs to search for
    2) input_file.txt = a file listing all the count files to search. Count
        files should be two columns, first column is spacer construct and second 
        column is count.
        
Output is a file of counts for each gene ID for each file in the input_file.txt

@author: kevinmyers
"""
import argparse
import re
import sys

def spacer_count(input_file):
    output = re.sub('_constructs_values.txt', '_spacer_counts.txt', input_file)
    id_list = []
    with open("spacer_ids.txt", 'r') as ids:
        for line in ids:
            line = line.rstrip()
            id_list.append(line)    
            
    gene_dict = {}
    with open(input_file, 'r') as f:
        for line in f:
            line = line.rstrip()
            construct = line.split('\t')[0]
            count = line.split('\t')[1]
            for each in id_list:
                if each in construct:
                    if each in gene_dict:
                        gene_dict[each].append(int(count))
                    else:
                        gene_dict[each]=[]
                        gene_dict[each].append(int(count))
                        
    out = {k: [sum(gene_dict[k])] for k in gene_dict.keys()}
    
    with open(output, 'w') as f:
        for key,val in out.items():
            f.write(f"{key}\t{val}\n")

def main():
   
    cmdparser = argparse.ArgumentParser(description="Calculate each spacer presence in files.",
                                        usage='%(prog)s -f <list of files to search> ' ,prog='individual_spacer_total_counts.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store',      dest='FILE' , help='File with files to process, one per line.', metavar='')
    cmdResults = vars(cmdparser.parse_args())
    
    search_files = []
    
    if cmdResults['FILE'] is not None:
        input_files = cmdResults['FILE']
    else:
        print("Please provide a files to search, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)
        
    with open(input_files, 'r') as f:
        for line in f:
            line2 = line.rstrip('\n')
            search_files.append(line2)
    
    for input_file in search_files:
        spacer_count ( input_file )
        
if __name__ == "__main__":
    main()