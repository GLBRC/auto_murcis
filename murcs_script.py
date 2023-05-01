#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
murcs_script.py

This script will count the number of times each spacer + repeat combination appears in FASTA files

Input is a list of the spacer and repeat combinations and list of FASTA files to search

Spacer/Repeat Combination File should be a two column, tab delimited file with the first tab
being the name of the Spacer/Repeat combindation and the second column the sequence to search (-t):
    
    ID                  Seq
    Gene_S_R_1          ATTAAGCCATGGCAGTGCAGACGATAGAGCACATAGCTAGCTATACGATAAAATCG
    Gene_S_R_2          TTCAAAACAGCATAGCTCTAAAACATTAAGCCCCAAAAAATAGATATGAGCCACAG
    
Fasta File should be a list with FASTA files to process in a single column, one per line (-f)

Path should be indicated as the location of the Rscript (-p)    

There are multiple outputs, each named by remove the .fasta and adding text:
    _search_results.txt = Initial search results file
    _gene_combinations_per_read.txt = list of genes with matches to each read
    _gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by gene name
    _gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries written to file (position and gene name)
    _gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted by position of match relative to the start of the read (spacer order)
    All_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible constructs and sub-constructs along with counts for each experiment
    temp_sorted_possible_combinations.txt = intermediate file with all possible construct and sub-construct combinations
    temp_sorted_possible_combinations_Dict.txt = intermediate file with dictionary of all possible construct and sub-construct combinations
    _pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with minHit of 5
    _pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample including all hits
    _pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample with minHit of 5
    
Requirements:
    Python modules:
        argparse
        SeqIO from BioPython
        Combinations from Itertools
        os
        re
        subprocess
        sys
        time
    R libraries:
        ggplot2
        reshape2
        circlize

Script must be run in the same directory as the FASTA files and Spacer/Repeat Combindation File
Requires Python3 and R

This has been tested on MacOS 12.6. It may require modification to run on other operating systems

@author: kmyers2@wisc.edu
"""
from Bio import SeqIO
from contextlib import redirect_stderr
from itertools import combinations
import argparse
import pyfastx
import pysam
import os
import re
import subprocess
import sys
import time

def search( gene_list, fasta ):
    """search

    Search for matches to every spacer+repeat combination and organize the
    output and indicate number of times each construct and sub-construct
    are found in the data.
    
    Parameters
    ----------
    gene_list : list
        text file with each possible spacer+repeat combination
    fasta : str 
        fasta file from PacBio CCS processing to search

    Returns
    -------
    _search_results.txt = Initial search results file
    _gene_combinations_per_read.txt = list of genes with matches to each read
    _gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by gene name
    _gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries written to file (position and gene name)
    _gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted by position of match relative to the start of the read (spacer order)

    """
    search_terms = {}
    out_searchResults = fasta.replace('.fasta', '_search_results.txt')
    out_geneCombo = fasta.replace('.fasta', '_gene_combinations_per_read.txt')
    out_geneCombo_sorted = fasta.replace('.fasta', '_gene_combinations_per_read_SortedByGeneName.txt')
    
    print(f"Working on {fasta} now…\n")

# Read in Spacer + Repeat list and write to dictionary

    with open(gene_list, 'r') as f:
        for _ in range(1):
            next(f)
        for line in f:
            (k, v) = line.split()
            search_terms[k]=v

# Using Biopython, read in sequencing FASTA file and write to dictionary
    
    seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(fasta, 'fasta')}

# Search for matches to Spacer + Repeat combinations in the FASTA file
# Report read name, gene (spacer) name, and position of the start of the 
# Spacer + Repeat matches realative to the start of the sequence.

    search_results = []
    for key,value in search_terms.items():
        for key2,value2 in seq_dict.items():
            if re.search(str(value), str(value2)):
                position = re.search(str(value), str(value2)).start()
                search_results.append(f"{key2}\t{key}\t{position}\n")
            else:
                continue
            #if value in value2:
            #   search_results.append(f"{key2}\t{key}\n")           
   
 # Write results to file for notes and for future use in the script   
   
    with open(out_searchResults, 'w') as f:
        f.write("Read_Name\tGenes_and_Positions\n")
        for each in search_results:
            f.write(each)

# Read in search file, parse it, and add to dictionary

    search_results_dict = {}
    with open(out_searchResults, 'r') as f:
        for _ in range(1):
            next(f)
        for line in f:
            geneCombo = line.split('\t')[1].rstrip()
            readName = line.split('\t')[0]
            geneID = geneCombo.split('_')[0]
            if readName in search_results_dict:
                search_results_dict[readName].append(f"{geneID}\t")
            else:
                search_results_dict[readName]=[]
                search_results_dict[readName].append(f"{geneID}\t")

# Combine all genes for each read in a list (remove dubplicates with set)
# write this to a file for notes
    
    for i in search_results_dict:
        search_results_dict[i] = list(set(search_results_dict[i]))
    with open(out_geneCombo, 'w') as f:
        for key, value in search_results_dict.items():
            f.write(f"{key}\t{value}\n")

# Go through the dictionary (read name and list of genes with matches) and 
# add to list, then edit to make it pretty. Write this to a file for notes.
 
    gene_combo_final = []
    for key, val in search_results_dict.items():
        gene_combo_final.append(f"{key}\t{val}\n")  
    gene_combo_final2 = []
    for line in gene_combo_final:
        line2 = re.sub("\]", "", line)
        line2 = re.sub("\[", "", line2)
        line2 = re.sub("['']", "", line2)
        line2 = line2.replace('\\t', '\t')
        line2 = line2.replace('\t, ', ', ')
        line2 = line2.replace('\t\n', '\n')
        gene_combo_final2.append(line2)    
    with open(out_geneCombo, "w") as f:
        f.write("Read_Name\tGenes\n")
        for each in gene_combo_final2:
            f.write(each)

# Sort the list based on gene name and write to a new list and then a file 
            
    sortedList = []        
    with open(out_geneCombo, 'r') as f:
        for _ in range(1):
            next(f)
        for line in f:
            seqName = line.split('\t')[0]
            genes = line.rstrip('\n').split('\t')[1]
            geneList = (genes.split(', '))
            geneList.sort()
            sortedList.append(f"{seqName}\t{geneList}\n")
    sortedList2 = []
    for line in sortedList:
        line2 = line.replace("[", "")
        line2 = line2.replace("]", "")
        line2 = line2.replace("'", "")
        sortedList2.append(line2)

    with open(out_geneCombo_sorted, 'w') as f:
        f.write("Read_Name\tGenes_Sorted_By_Gene_Name\n")
        for each in sortedList2:
            f.write(each)

# Now do the same search, but this time include the position first to sort based
# on position from the start of the read.
    
    out_dictionary_to_file = fasta.replace('.fasta', '_gene_combinations_per_read_dictionary_out.txt')
    out_geneCombo_sorted_final = fasta.replace('.fasta', '_gene_combinations_per_read_SortedBySpacerPosition.txt')
    
# Search for matches to Spacer + Repeat combinations in the FASTA file
# Report read name, gene (spacer) name, and position of the start of the 
# Spacer + Repeat matches realative to the start of the sequence.
        
    search_results_with_position = {}        
    for search_term, search_seq in search_terms.items():
        for read_name,read_seq in seq_dict.items():        
            if re.search(str(search_seq), str(read_seq)):
                position = re.search(str(search_seq), str(read_seq)).start()
                gene_ID = search_term.split('_')[0]
                position_geneID_tuple = (str(position), str(gene_ID))
                position_geneID = ','.join(position_geneID_tuple)
                if (read_name) in search_results_with_position:
                    if (gene_ID) in search_results_with_position[read_name]:
                        continue
                    else:
                        search_results_with_position[read_name][gene_ID] = []
                        search_results_with_position[read_name][gene_ID].append(position_geneID)
                else:
                    search_results_with_position[read_name]={}
                    if (gene_ID) in search_results_with_position[read_name]:
                        continue
                    else:
                        search_results_with_position[read_name][gene_ID] = []
                        search_results_with_position[read_name][gene_ID].append(position_geneID)
            else:
                continue

# Write the dictionary to a file for notes

    with open(out_dictionary_to_file, 'w') as f:
        for key, val in search_results_with_position.items():
            f.write(f"{key}\t{val}\n")

# Read in dictionary file and write to a list
            
    pos_gene_list = []
    with open(out_dictionary_to_file, 'r') as f:
        for line in f:
            seqName = line.split('\t')[0]
            genes = line.rstrip('\n').split('\t')[1]
            pos_gene_list.append(f"{seqName}\t{genes}\n")

# Clean up the list to make it pretty and sortable by postion
    
    Pos_gene_list_edited = []
    for line in pos_gene_list:
        line2 = line.replace("{", "")
        line2 = line2.replace("}", "")
        line2 = line2.replace("'", "")
        line2 = re.sub("lpg[0-9]*: ", "", line2)
        line2 = line2.replace("[", "")
        line2 = line2.replace("]", "")
        Pos_gene_list_edited.append(line2)

# Sort the list by position for each spacer match
    
    pos_gene_sorted_list = []
    pos_gene_sorted_list_edited = []
    for line in Pos_gene_list_edited:
        seqName = line.split('\t')[0]
        genes = line.rstrip('\n').split('\t')[1]
        geneList = (genes.split(', '))
        geneList.sort()
        pos_gene_sorted_list.append(f"{seqName}\t{geneList}\n")

# Clean up the list again, getting rid of position and leave just a list of
# gene names, sorted by the position of the match relative to the start
# of the read. Write this to a file.
        
    for line in pos_gene_sorted_list:
        line2 = re.sub("'[0-9]*,l", "l", line)
        line2 = line2.replace("]", "")
        line2 = line2.replace("[", "")
        line2 = line2.replace("'", "")
        pos_gene_sorted_list_edited.append(line2)
        
    with open(out_geneCombo_sorted_final, 'w') as f:
        f.write("Read_Name\tGene_Spacer_Sorted_By_Position_in_Read\n")
        for each in pos_gene_sorted_list_edited:
            f.write(each)
            
# Take the sorted by gene name results, import and for all constructs >2,
# determine all possible combinations (no replacement) and then count how often
# each construct and sub-construct has a match. Write to file

# Import all constructs and write to list

    gene_groups = []
    with open(out_geneCombo_sorted, 'r') as f:
        for line in f:
            gene_groups.append(line.rstrip('\n').split('\t')[1])

# Select only those constructs with >2 spacers and write to new list
            
    greater_than_2plex = []
    for each in gene_groups:
        if len(each) > 16:
            greater_than_2plex.append(each)
    greater_than_2plex = list(set(greater_than_2plex))

# For each construct >2 spacers, write to dictionary key and determine all the 
# possible combinations (no replacement) and wite as value in decending order
# based on construct length. Then write the dictionary to a list. 

    combo_Plex_Dict = {}
    for each in greater_than_2plex:
        combo_Plex_Dict[each]=[]
        eachList = each.split(', ')
        length = len(eachList)
        plexRange = range(length-1, 0, -1)
        for i in plexRange:
            combo_Plex_Dict[each].append(list(combinations(eachList, i)))

    combo_Plex = []
    for key, value in combo_Plex_Dict.items():
        combo_Plex.append(key)
        for each in value:
            combo_Plex.append(each)
            
# Write the intermediate result to a file (for records)

    out_spacer_constructs_sub_first = fasta.replace('.fasta', '_all_spacer_constructs_and_subconstructs.txt')
    
    with open(out_spacer_constructs_sub_first, 'w') as f:
        for key, value in combo_Plex_Dict.items():
            f.write(f"\n{key}\n{value}\n")

# Read in the constructs and sub-constructs file, clean it up and write to 
# new list and new file (for records)

    cleaned_construct_subconstruct_list = []        
    with open(out_spacer_constructs_sub_first, 'r') as f:
        for line in f:
            line2 = line.replace(", (", "\n")
            line2 = line2.replace(", [(", "\n")
            line2 = line2.replace("'", "")
            line2 = line2.replace("[", "")
            line2 = line2.replace("(", "")
            #line2 = line2.replace(", ", ",")
            line2 = line2.replace("]", "")
            line2 = line2.replace(")", "")
            line2 = line2.replace(",\n", "\n")
            cleaned_construct_subconstruct_list.append(line2)
    
    out_spacer_constructs_sub_cleaned = fasta.replace('.fasta', '_all_spacer_constructs_and_subconstructs_organized.txt')        
    
    with open(out_spacer_constructs_sub_cleaned, 'w') as f:
        for each in cleaned_construct_subconstruct_list:
            f.write(each)

# Read in cleaned and orgainzed construct and sub-construct file and write to
# list. Remove the blank line first.
            
    plex_to_search = []
    with open(out_spacer_constructs_sub_cleaned, 'r') as f:
        for line in f:
            if line.startswith("l"):
                line2 = line.rstrip('\n')
                plex_to_search.append(line2)
            else:
                continue

# For each construct and sub-construct, search for an exact match in the original
# gene_groups list and count the number of matches using len(matching)
# Then write it all to a final file.
            
    matching_list = []
    for each in plex_to_search:
        matching = [s for s in gene_groups if each == s]
        number = len(matching)
        matching_list.append(f"{each}\t{number}\n")
    
    out_spacer_constructs_sub_final = fasta.replace('.fasta', '_all_spacer_constructs_and_subconstructs_with_number_of_matches.txt')     
        
    with open(out_spacer_constructs_sub_final, 'w') as f:
        f.write("Spacer_Construct\tNumber_of_Hits\n")
        for each in matching_list:
            f.write(each)

# Determine the counts for each pair of spacers
# Start by reading in the original spacer file, parsing it and writing the
# combinations (without replacement) to a new list and sorting

    spacer_id = []
    with open(gene_list, 'r') as f:
        f.readline()
        for line in f:
            spacerID = line.split('\t')[0]
            spacer_name = spacerID.split('_')[0]
            spacer_id.append(spacer_name)
    spacer_id = list(set(spacer_id))
    
    spacer_id.sort()
    
    spacer_id_combo = list(combinations(spacer_id, 2))
    
    spacer_id_combo_string = []
    for each in spacer_id_combo:
        spacer_id_combo_string.append(str(each))

# Clean up the spacer combo list
    
    updated_spacer_id_combo = []
    
    for line in spacer_id_combo_string:
        line2 = line.replace("('", "")
        line2 = line2.replace("', '", ", ")
        line2 = line2.replace("')", "")
        updated_spacer_id_combo.append(line2)

# Open the sorted by spacer name with value file and write the constructs to
# a new list for searching.
    
    in_file_sorted = []
    with open(out_geneCombo_sorted, 'r') as f:
        f.readline()
        for line in f:
            grow = line.rstrip('\n').split('\t')[1]
            if len(grow) > 9:
                spacer_name = grow.split(', ')
                in_file_sorted.append(list(combinations(spacer_name, 2)))
                
    in_file_sorted_to_search = []
    for each in in_file_sorted:
        for tupID in each:
            strID = ', '.join(tupID)
            in_file_sorted_to_search.append(strID)

# Find matches to paired spacers and write to a file with the vaule
        
    pair_matching_list = []
    
    for each in updated_spacer_id_combo:
        matching = [s for s in in_file_sorted_to_search if each == s]
        number = len(matching)
        pair_matching_list.append(f"{each}\t{number}\n")
        
    out_pairwise_without_replacement_value = fasta.replace('.fasta', '_pairwise_spacer_combinations_withoutReplacement_value_for_Chord_Diagrams.txt')
    
    with open(out_pairwise_without_replacement_value, 'w') as f:
        f.write("Spacer_1, Spacer_2\tCount\n")
        for each in pair_matching_list:
            f.write(each)
    
    out_pairwise_without_replacement_value_forPlotting_list = []
    
    out_pairwise_without_replacement_value_forPlotting = fasta.replace('.fasta', '_pairwise_spacer_combinations_withoutReplacement_value_for_Chord_Diagrams_forPlotting.txt')
    
    with open(out_pairwise_without_replacement_value, 'r') as f:
        for line in f:
            line2 = line.replace(", ", "\t")
            out_pairwise_without_replacement_value_forPlotting_list.append(line2)
            
    with open(out_pairwise_without_replacement_value_forPlotting, 'w') as f:
        for each in out_pairwise_without_replacement_value_forPlotting_list:
            f.write(each)
            
# Take the list of spacer pairs and write a new file that contains both 
# combintations of the pair (gene 1-gene2 and gene2-gene1) and the value
# as well as the same spacer pair (gene1 - gene1) and the artifical value of
# 10,000. This is for the correlation plots.
            
    all_pairs_list = []
    with open(out_pairwise_without_replacement_value, 'r') as f:
        f.readline()
        for line in f:
            construct = line.split('\t')[0]
            gene1 = construct.rstrip().split(', ')[0]
            gene2 = construct.rstrip().split(', ')[1]
            value = line.rstrip('\n').split('\t')[1]
            all_pairs_list.append(f"{gene1}, {gene2}\t{value}\n")
            all_pairs_list.append(f"{gene2}, {gene1}\t{value}\n")
    for each in spacer_id:
        all_pairs_list.append(f"{each}, {each}\t10000\n")
        
    all_pairs_list.sort()
    
    out_pairwise_WITH_replacement_value = fasta.replace('.fasta', '_pairwise_spacer_combinations_withReplacement_value_for_correlation_plots.txt')
    
    with open(out_pairwise_WITH_replacement_value, 'w') as f:
        f.write("Spacer_1, Spacer_2\tCount\n")
        for each in all_pairs_list:
            f.write(each)
            
    out_pairwise_WITH_replacement_value_forPlotting_list = []
    
    out_pairwise_WITH_replacement_value_forPlotting = fasta.replace('.fasta', '_pairwise_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt')
    
    with open(out_pairwise_WITH_replacement_value, 'r') as f:
        for line in f:
            line2 = line.replace(", ", "\t")
            out_pairwise_WITH_replacement_value_forPlotting_list.append(line2)
            
    with open(out_pairwise_WITH_replacement_value_forPlotting, 'w') as f:
        for each in out_pairwise_WITH_replacement_value_forPlotting_list:
            f.write(each)
    

def combineCountFiles( cwd ):
    """combineCountFiles

    Organize every possible construct and sub-construct from all experiments and
    record the hits for each in a single file.
    
    Parameters
    ----------
    cwd : str
        Current working directory
    
    Returns
    -------
    All_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible constructs and sub-constructs along with counts for each experiment
    temp_sorted_possible_combinations.txt = intermediate file with all possible construct and sub-construct combinations
    temp_sorted_possible_combinations_Dict.txt = intermediate file with dictionary of all possible construct and sub-construct combinations

    """
    all_possible_constructs_greaterThan_2 = []

    individual_count_files = [ fn for fn in os.listdir(cwd) if fn.endswith("constructs_and_subconstructs_with_number_of_matches.txt") ]
    individual_count_files.sort()

    for each in individual_count_files:
        with open(each, 'r') as f:
            for _ in range(1):
                next(f)
            for line in f:
                construct = line.rstrip('\n').split('\t')[0]
                if len(construct) > 16:
                    all_possible_constructs_greaterThan_2.append(construct)
                    
    sorted_construct_list = []

    sorted_construct_list = list(set(all_possible_constructs_greaterThan_2))
    sorted_construct_list = sorted(sorted_construct_list, key = len, reverse = True)
            
    sorted_possible_combinations_Dict = {}
    for each in sorted_construct_list:
        sorted_possible_combinations_Dict[each]=[]
        eachList = each.split(', ')
        length = len(eachList)
        plexRange = range(length-1, 0, -1)
        for i in plexRange:
            sorted_possible_combinations_Dict[each].append(list(combinations(eachList, i)))
            
    combo_Plex = []
    for key, value in sorted_possible_combinations_Dict.items():
        combo_Plex.append(key)
        for each in value:
            combo_Plex.append(each)
            
    # Write the intermediate result to a file (for records)

    with open('temp_sorted_possible_combinations_Dict.txt', 'w') as f:
        for key, value in sorted_possible_combinations_Dict.items():
            f.write(f"\n{key}\n{value}\n")

    # Read in the constructs and sub-constructs file, clean it up and write to 
    # new list and new file (for records)

    cleaned_sorted_possible_combo_list = []        
    with open('temp_sorted_possible_combinations_Dict.txt', 'r') as f:
        for line in f:
            line2 = line.replace(", (", "\n")
            line2 = line2.replace(", [(", "\n")
            line2 = line2.replace("'", "")
            line2 = line2.replace("[", "")
            line2 = line2.replace("(", "")
            #line2 = line2.replace(", ", ",")
            line2 = line2.replace("]", "")
            line2 = line2.replace(")", "")
            line2 = line2.replace(",\n", "\n")
            cleaned_sorted_possible_combo_list.append(line2)

    with open('temp_sorted_possible_combinations.txt', 'w') as f:
        for each in cleaned_sorted_possible_combo_list:
            f.write(each)
    
    
    count_files = [ fn for fn in os.listdir(cwd) if fn.endswith("spacer_constructs_and_subconstructs_with_number_of_matches.txt")]
    count_files.sort()
    combined_files_Dict = {}
    count = 0
    maxLen = 0
    for c, each in enumerate(count_files, 1):
        with open(each, 'r') as f:
            temp_construct_list = []
            f.readline()
            for line in f:
                construct_count = line.rstrip('\n')
                temp_construct_list.append(construct_count)
                set_construct_list = set(temp_construct_list)
            for each in set_construct_list:
                construct = each.split('\t')[0]
                NumHits = each.rstrip('\n').split('\t')[1]
                if construct not in combined_files_Dict:
                    combined_files_Dict[construct]=[] 
                    if count > 0:
                        for i in range(count):
                            combined_files_Dict[construct].append('0')
                        combined_files_Dict[construct].append(NumHits)
                    else:
                        combined_files_Dict[construct].append(NumHits) 
                else:
                    combined_files_Dict[construct].append(NumHits)
            count += 1        
        maxLen = max({len(x) for x in combined_files_Dict.values()})
        
        for k,v in combined_files_Dict.items():
            if len(v) < maxLen:
                v.append('0')
                
    sample_names = []
    for each in count_files:
        sampleName = each.split('.ccs')[0]
        sample_names.append(sampleName)
        
    sampleNameHeader = '\t'.join(sample_names)
    
    possible_construct_list = []
    with open('temp_sorted_possible_combinations.txt', 'r') as f:
        f.readline()
        for line in f:
            line2 = line.rstrip('\n')
            possible_construct_list.append(line2)

    new_list = []                
    for each in possible_construct_list:
        for key, val in combined_files_Dict.items():
            if each == key:
                new_list.append(f"{each}\t{val}\n")
                
    final_list = []
    
    for line in new_list:
        line2 = line.replace("[", "")
        line2 = line2.replace("]", "")
        line2 = line2.replace(", '", "\t")
        line2 = line2.replace("'", "")
        final_list.append(line2)
        
    with open('All_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt', 'w') as f:
        f.write(f"Spacer_Consturct\t{sampleNameHeader}\n")
        for each in final_list:
            f.write(f"{each}")

def chord_correlation_plots( cwd, path_to_R ):
    """chord_correlation_plots

    Run the Rscript to plot the Chord and Correlation plots
    
    command = Rscript /Users/kevinmyers/scripts/nicole_nih_scripts/chord_correlation_plot_script.R ./
    
    Output:
        Chord diagrams for all samples (minHit5 only)
        Correlation charts for all samples (all and minHit5 only)
    """    
    program = path_to_R + '/chord_correlation_plot_script.R'
    cmd = [ 'Rscript', program , cwd ]
    subprocess.run( cmd )  

def cleanUp( cwd ):
    """cleanUp

    Organize the working directory by moving all intermediate files to a new 
    directory (other_files) and the initial FASTA files to a new directory
    (fasta_files)

    Parameters
    ----------
    cwd : str 
        Current Working Directory

    """
    cwd = cwd + "/"
    # move intermediate files to another folder
    os.mkdir( "other_files" )
    bamDir = cwd + "/other_files/"
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("match_search_results.txt") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_per_read.txt") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_dictionary_out.txt") ] 
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("constructs_and_subconstructs.txt") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_constructs_and_subconstructs_organized.txt") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.startswith("temp") ]
    # organize the original fasta files
    os.mkdir( "fasta_files" )
    bamDir = cwd + "/fasta_files/"
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".fasta") ]

def makeFasta(bamList):
    """makeFasta

    Using list of bam files create a new fasta file for each bam.

    Parameters
    ----------
    bamList : list
        List containing bam files to process.
    Returns
    -------
    fastaLst : list
        List of newly created fasta files.
    
    """
    fastaLst = []           # list of newly created fasta files to return

    with open('bamToFasta.log', 'w') as log:                  # log any errors
    # create a fasta file for each bam file using pysam
        for bam in bamList:
            fastaName = re.sub('.bam', '.fasta',bam)             # create new fasta name
            cmd = ['samtools', 'fasta', '-0', fastaName, bam] # setup samtools command
            # run command
            output = subprocess.Popen(cmd, stderr=subprocess.PIPE).communicate()
            # log any problems
            result = output[1].decode('utf-8')
            log.write(result)
            log.write("\n")
    log.close()

    return fastaLst 

def countLines(fasta):
    """countLines
    
    Count lines in fasta file and divide by 2 to return the number of sequences.
    Utilizes Pysam to get sequence lengths and number of reads.

    Parameters
    ----------
    fasta : str
        fasta file name
    Returns
    -------
    count : int
        Number of sequences in fasta file
    """
    fsa = pysam.FastaFile(fasta)
    seqLens  = fsa.lengths
    numReads = len(fsa.lengths)
    return seqLens, numReads    

def processFasta(FastaLst):
    """processFasta

    Collect read counts, read lengths for all fasta files and write to file.

    Parameters
    ----------
    FastaLst : list
        List containing the fasta files to process (previously converted from input bams)    
    """    
    pass   


        
    
def main():
   
    cmdparser = argparse.ArgumentParser(description="Count spacers + repeat in FASTA files for NIH Project and produce organized files for further analysis along with different plots.",
                                        usage='%(prog)s -f <list of FASTAs to process> -t <spacer and repeat combinations> -p <path to GitHub Repo and Rscript> [optional arguments: -d]', prog='count_spacers_NIH.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE' , help='File with all the bam files to process, one per line.', metavar='')
    cmdparser.add_argument('-t', '--targets',  action='store', dest='TARGETS' , help='spacer + repeat targets for each gene', metavar='')
    cmdparser.add_argument('-p', '--path', action='store', dest='PATH', help='Path to location of GitHub Repo and Python and R Scripts', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL' , help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
    
    cwd   = os.getcwd()
    
    start = time.time() #start timer to see how long the program takes to run
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

    if cmdResults['DETAIL']:
        print("\nmurcs_script.py")
        print("\nPurpose: Count spacers + repeat in FASTA files for NIH Project and produce organized files for further analysis along with different plots.")
        print("\nInput: A text file with each FASTA to process and a list of all spacer+repeat sequences and names.")
        print("\nPlease use a dedicated directory for running this pipeline.")
        print("Generate the directory and copy the FASTA and spacer+repeat files to the new directory.")
        print("Produce FASTA input file by running ls *fasta > fasta_input.txt")
        print("Required parameters:")
        print("\t1) -f fasta_input.txt")
        print("\t2) -t spacer_repeat_seq.txt\n")
        print("\t3) -p path_to_scripts\n")
        print("Optional arguments:")
        print("\t-d:  print a detailed description of the program.\n\n")
        print("Intermediate files are moved to the 'other_files' directory when finished")
        print("\nOriginal FASTA files are moved to the 'fasta_files' directory when finished")
        print("\nThere are multiple outputs, each named by remove the .fasta and adding text:")
        print("\t_search_results.txt = Initial search results file")
        print("\t_gene_combinations_per_read.txt = list of genes with matches to each read")
        print("\t_gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by gene name")
        print("\t_gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries written to file (position and gene name")
        print("\t_gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted by position of match relative to the start of the read (spacer order)")
        print("Three output files are for all experiments combined:")
        print("\tAll_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible constructs and sub-constructs along with counts for each experiment")
        print("\temp_sorted_possible_combinations.txt = intermediate file with all possible construct and sub-construct combinations")
        print("\ttemp_sorted_possible_combinations_Dict.txt = intermediate file with dictionary of all possible construct and sub-construct combinations")
        print("\t_pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with minHit of 5")
        print("\t_pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample including all hits")
        print("\t_pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample with minHit of 5")
        print("\n")
        print("See Kevin Myers (kmyers2@wisc.edu) with any questions.\n\n")
        sys.exit(1)
    
    BAM_files   = []   # hold list of initial input bam files to process
    
    if cmdResults['FILE'] is not None:
        bamfile = cmdResults['FILE']
        with open(bamfile, 'r') as f:
            for bam in f:
                BAM_files.append(bam.rstrip())
    else:
        print("Please provide a file with GFF names, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)

    # create fasta files from bam file, store in list
    FASTA_files = makeFasta(BAM_files)#################################
    ## TESTING ONLY
    #Fasta_files = ['Amoeba_output_A2.ccs.fasta']

    # retrieve spacer and repeat file
    if cmdResults['TARGETS'] is not None:
        gene_list = cmdResults['TARGETS']
    else:
        print("Please provide a file spacer + repeat sequences.\n")
        cmdparser.print_help()
        sys.exit(1)

    '''    
    if cmdResults['PATH'] is not None:
        path = cmdResults['PATH']
    else:
        print("Please provide a path to the Python and R scripts (GitHub Repo).\n")
        cmdparser.print_help()
        sys.exit(1)
    '''
    # report number of fasta files to process
    number_of_files = len(Fasta_files)
    print(f"There are {number_of_files} Fasta files to process.\n")    
    '''
    # create an output file containing all the read lengths of the original input files 
    readStats    = {}          # store a list of all lengths for all samples    
    orig_readLen = []          # 
    for fsa in Fasta_files:
        seqLen, totalReads = countLines(fsa)
        if fastq not in readStats:
            readStats[fastq] = seqLen
        with open(re.sub('.fasta', '_read_lengths.txt',fsa)) as out:
            for rd in totalReads:
                out.write(f'{rd}\n')
        out.close()
    '''    
    '''
    for fasta in FASTA_files:
        search( gene_list, fasta )
    
    print("Combining all the experimental count files together…\n")
    
    combineCountFiles( cwd )
    
    print("Plotting the Chord and Correlation plots…\n")
    
    if path.endswith('/'):
        path_to_R = path[:-1]  #remove last / from path
    else:
        path_to_R = path
        
    chord_correlation_plots( cwd, path_to_R )
    
    print("Let's clean up and get out of here!\n")
    
    cleanUp( cwd )
    
    # end timer and do math and report how long the script took to run
    end = time.time()
    total_time = round(end - start, 2)
    total_time_min = round(total_time/60, 2)
    total_time_hours = round(total_time/60/60, 2)
    print(f"\nIt took {total_time_hours} hours ({total_time_min} minutes) to process the {number_of_files} FASTA files.\n")
    print("Please email Kevin Myers (kmyers2@wisc.edu) with any questions.\n")
    '''
        
if __name__ == "__main__":
    main()