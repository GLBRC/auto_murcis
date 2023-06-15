#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
murcs_script.py

This script will count the number of times each spacer + repeat combination appears in FASTA files

Input
-----

Spacer/repeat sequence file should be two columns, tab delimited file with the first column
being the name of the Spacer/Repeat and the second column the sequence to search:
    
    name        sequence
    repeat      GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC
    lpg0404a	GCAAGTTACGGATACTGTCTACAG
    lpg0404b	GGGTGCCTTGTTCTTCCTTATCCC
    
Bam list file should be a bam file names to process in a single column, one per line (-f)  

Output
------

There are multiple outputs, each named by remove the .fasta and adding text:
    _search_results.txt = Initial search results file
    _gene_combinations_per_read.txt = list of genes with matches to each read
    _gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by 
    gene name
    _gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries
      written to file (position and gene name)
    _gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted
      by position of match relative to the start of the read (spacer order)
    All_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible 
      constructs and sub-constructs along with counts for each experiment
    temp_sorted_possible_combinations.txt = intermediate file with all possible
      construct and sub-construct combinations
    temp_sorted_possible_combinations_Dict.txt = intermediate file with dictionary
      of all possible construct and sub-construct combinations
    _pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with 
      minHit of 5
    _pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample 
      including all hits
    _pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample 
      with minHit of 5
    
Requirements
------------
    Python 3
    Python modules:
        BioPython
    R
    R libraries:
        ggplot2
        reshape2
        circlize

Script must be run in the same directory as the bam and Spacer/Repeat Combindation files.

Script tested on MacOS 12.6 and CentOS Linux 7(Core). 
It may require modification to run on other operating systems

@author: kmyers2@wisc.edu
"""
from Bio.Seq import Seq
from Bio import SeqIO
from datetime import date
from functools import reduce
import argparse
import glob 
import itertools                                  # for zip_longest and combinations
import logging 
import multiprocessing as mp
import pandas as pd                                # use data frames for merging counts
import pyfastx
import pysam
import os
import pwd
import re
import subprocess
import sys
import time

# set path to script home
dirPath = (os.path.dirname(os.path.abspath(__file__))) + '/'

def countGenes():
    """countGenes

    Count genes (spacer/Repeat hits) for each sample.

    Returns
    -------
    pandas data frame    
    """
    cwd = os.getcwd()
    dfLst = []

    # get a list of files to process, by globbing *-subset_gene_combinations_per_read.txt files
    for file in glob.glob(cwd + '/*-gene_combinations_per_read.txt'):
        cnt = {}   # key = spacer/Repeat value = count
        sample = os.path.basename(re.sub('-gene_combinations_per_read.txt', '', file))
        with open(file, 'r') as f:
            f.readline()   # skip header
            for ln in f:
                read, genes = ln.rstrip().split('\t')
                if genes not in cnt:
                    cnt[genes] = 1
                else:
                    cnt[genes] += 1

        df = pd.DataFrame.from_dict(cnt, orient='index', columns=[sample])
        dfLst.append(df)
        # write individual sample counts to file
        df.to_csv( sample + '-spacerRepeat-CNTs.txt', sep='\t', index_label='spacer/Repeat', na_rep=0)

    data_merge = reduce(lambda left,right: pd.merge(left,right, how="outer", left_index=True, right_index=True),dfLst)
    
    return data_merge

def ngrams(seq, n=60):
    """ngrams
    
    Create a list of ngrams using a DNA sequence string.

    Parameters
    ----------
    seq : str
        DNA sequence 
    n : int
        ngram size, defaults to 60.

    Returns
    -------
    seqNgrams : list
        list of all ngrams in seq.  

    ngram recipe from: https://bergvca.github.io/2017/10/14/super-fast-string-matching.html
    """
    ngrams = zip(*[seq[i:] for i in range(n)])
    return [''.join(ngram) for ngram in ngrams]

def search(geneTags, fasta, sizeSprRpt):
    """search

    Search for matches to every spacer+repeat combination and organize the
    output and indicate number of times each construct and sub-construct
    are found in the data.
    
    Parameters
    ----------
    geneTags : list
        list with each possible spacer+repeat combination
    fasta : str 
        fasta file from PacBio CCS processing to search
    sizeSprRpt : int
        Expected length of space repeat.  
    Returns
    -------
    _search_results.txt = Initial search results file
    _gene_combinations_per_read.txt = list of genes with matches to each read
    _gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by gene name
    _gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries written to file (position and gene name)
    _gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted by position of match relative to the start of the read (spacer order)

    """
    search_terms = {}
    out_searchResults = fasta.replace('.fasta', '-search_results.txt')
    out_geneCombo     = fasta.replace('.fasta', '-gene_combinations_per_read.txt')
    out_geneOrder     = fasta.replace('.fasta', '-gene-order_byRead.txt')
    
    print(f"Working on {fasta} now…\n")
    # open fasta file to process and output files 
    with open(fasta, 'r') as f, open(out_searchResults, 'w') as sr, open(out_geneCombo,'w') as cout, open(out_geneOrder, 'w') as gout:
        sr.write('Read_name\tGenes\tPositions\n')
        cout.write('Read_name\tGenes\n')
        gout.write('Read_name\tGene_order\n')
        geneCombos = {}                      # collate all genes for a particular read
        geneOrder = {}                       # key= read_name value= list of spacerRepeats
        for header, read in itertools.zip_longest(*[f]*2):
            header = re.sub('>', '', header.rstrip())                 # read name 
            read   = read.rstrip()
            # search for ngrams 
            for pos, ngram in enumerate(ngrams(read, sizeSprRpt)):
                if ngram in geneTags:
                    outLine = f'{header}\t{geneTags[ngram]}\t{pos}\n'
                    sr.write(outLine)
                    # gather reads with gene names
                    if header not in geneCombos:
                        geneCombos[header] = set()
                        geneCombos[header].add(geneTags[ngram].split('_')[0])
                    else:
                        geneCombos[header].add(geneTags[ngram].split('_')[0])
                    # gather reads and full spacer/Repeat name
                    if header not in geneOrder:
                        geneOrder[header] = []
                        geneOrder[header].append(geneTags[ngram].split('_')[0])
                    else:
                        geneOrder[header].append(geneTags[ngram].split('_')[0])
                   
        # write gene combinations by read to file
        for readHdr in geneCombos.keys():
            cout.write(f'{readHdr}\t{",".join(list(geneCombos[readHdr]))}\n')
        # write gene order by read information to file
        for readHdr in geneOrder.keys():
            # remove gene name duplicates while maintaining order, using a dictionary
            geneLst = list(dict.fromkeys(geneOrder[readHdr]))  # get an ordered uniq list of gene names
            gout.write(f'{readHdr}\t{",".join(geneLst)}\n')   

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

def chord_correlation_plots(path_to_R):
    """chord_correlation_plots

    Run the Rscript to plot the Chord and Correlation plots
    
    command = Rscript /Users/kevinmyers/scripts/nicole_nih_scripts/chord_correlation_plot_script.R ./
    
    Output:
        Chord diagrams for all samples (minHit5 only)
        Correlation charts for all samples (all and minHit5 only)
    """    
    cwd = os.getcwd()
    program = path_to_R + '/chord_correlation_plot_script.R'
    cmd = ['Rscript', program, cwd]
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
            fastaName = re.sub('.bam', '.fasta',bam)          # create new fasta name
            fastaLst.append(fastaName)
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

def countRepeats(fsaFile, repeat):
    """countRepeats

    Determine the number of reads with repeat sequences, testing fwd & rvs orientations.
    Separate reads into new file for further processing & collect read lengths.

    Parameters
    ----------
    FastaLst : list
        List containing the fasta files to process (previously converted from input bams)    
    repeats : list
        Repeat list, used for searching.
    """    
    #for fsaFile in FastaLst:
    repeatsFile = re.sub('.fasta', '-repeat_match.fasta', fsaFile)
    readLenFile = re.sub('.fasta', '-repeats_match_read_lengths.txt', fsaFile)
    dat = pyfastx.Fasta(fsaFile,full_index=True)    # create a fasta object
    # open file to write the reads which contain the repeats, (filter out reads w/o repeat i.e. junk)
    with open(repeatsFile, 'w') as out, open(readLenFile, 'w') as lout:
        for i in range(dat.count(1)):           
            for rpt in repeat:
                if dat[i].search(rpt):
                    out.write(f'>{dat[i].name}\n')
                    out.write(f'{dat[i].seq}\n')
                    lout.write(f'{len(dat[i].seq)}\n')
    out.close()                 
    
def geneSpacerCombinations(spcrRpt):
    """geneSpacerCombinations

    Generate the 8 gene spacer combinations required to search the fasta files.
    These include Spacer+repeat, repeat+spacer, SpacerRC+repeat, repeat+SpacerRC,
    Spacer+repeatRC, repeatRC+spacer, SpacerRC+repeatRC, repeatRC+SpacerRC
    where RC indicates the reverse complement.

    Parameters
    ----------
    spcrRpt : str
        Spacer/repeat sequence file provided by user.
    Returns
    -------
    spacerRepeats : dictionary
        key = spacer/repeat value = gene name with spacer repeat order
    """
    spacerRepeats = {}
    with open(spcrRpt, 'r') as f, open('gene_spacer_combinations.txt', 'w') as out:
        f.readline()                                 # skip header
        rpt = f.readline().rstrip().split('\t')[1]   # grab repeat, 1st row
        rptRC = str(Seq(rpt).reverse_complement())
        
        for line in f:
            gene, spacer = line.rstrip().split('\t')
            spacerRC = str(Seq(spacer).reverse_complement())    
            
            # generate spacer repeat combinations
            #gene_Spacer+repeat
            name   =  gene + '_spacer+repeat' 
            target = spacer + rpt
            combo  = name + target
            out.write(combo + '\n')
            if target not in spacerRepeats:
                spacerRepeats[target] = name
            
            #gene_repeat+spacer  
            name   = gene + '_repeat+spacer'
            target = rpt + spacer
            combo  = name + target
            out.write(combo + '\n')
            if target not in spacerRepeats:
                spacerRepeats[target] = name            
            
            #gene_SpacerRC+repeat
            name   = gene + '_spacerRC+repeat'
            target = spacerRC + rpt
            combo  = name + target
            out.write(combo + '\n')
            if target not in spacerRepeats:
                spacerRepeats[target] = name            
            
            #gene_repeat+SpacerRC     
            name   = gene + '_repeat+spacerRC'
            target = rpt + spacerRC
            combo  =  name + target
            out.write(combo + '\n')            
            if target not in spacerRepeats:
                spacerRepeats[target] = name
            
            #gene_Spacer+repeatRC 
            name   = gene + '_spacer+repeatRC'    
            target = spacer + rptRC
            combo =  name + target
            out.write(combo + '\n')            
            if target not in spacerRepeats:
                spacerRepeats[target] = name
            
            #gene_repeatRC+spacer 
            name   = gene + '_repeatRC+spacer'   
            target = rptRC + spacer
            combo  = name + target
            out.write(combo + '\n')            
            if target not in spacerRepeats:
                spacerRepeats[target] = name
            
            #gene_SpacerRC+repeatRC  
            name   = gene + '_spacerRC+repeatRC'
            target = spacerRC + rptRC
            combo  = name + target
            out.write(combo + '\n')
            if target not in spacerRepeats:
                spacerRepeats[target] = name
            
            #gene_repeatRC+SpacerRC  
            name   = gene + '_repeatRC+spacerRC'
            target = rptRC + spacerRC
            combo  = name + target
            out.write(combo + '\n')            
            if target not in spacerRepeats:
                spacerRepeats[target] = name
    f.close()
    out.close()
    return spacerRepeats

def getRepeat(gene_list):
    """getRepeat

    Extract the repeat sequence from the spacer sequence and repeat file.

    searching for a line that looks like:

    repeat	GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC
    
    Parameters
    ----------
    gene_list : str
        File containg the repeat sequence.
    
    """
    fwd = None
    rvs = None 
    with open(gene_list, 'r') as f:
        f.readline()              # skip header
        for ln in f:
            if ln.startswith('repeat'):
                fwd = ln.rstrip().split('\t')[1]
                rvs = str(Seq(fwd).reverse_complement())
                break
    f.close()
    repeat = [fwd, rvs]
    return repeat        
    
def makePairwiseCnt():
    """makePairsizeCnt

    Generate 2 gene count tables 1) with all unidirectional pairs and 2) with each pair 
    represented twice, i.e.  gene1,gene2 and another pair gene2,gene1 with the same count.

    file expected to look like: 

    spacer/Repeat   Macro_output_B2.ccs
    lpg0026b,lpg0049a,lpg0228a,lpg1917b,lpg1917a    38
    lpg0049a,lpg0026b,lpg0228a      422
    lpg0281b,lpg0228a,lpg0228b,lpg0281a     580
    lpg1658a        200
    
    """
    for file in glob.glob( os.getcwd() + '/*-spacerRepeat-CNTs.txt'):
        fname = re.sub('-spacerRepeat-CNTs.txt', '', file) +  '_spacer_combinations_withoutReplacement_value_for_Chord_Diagrams_forPlotting.txt'
        fname2 = re.sub('-spacerRepeat-CNTs.txt', '', file) + '_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt'

        with open(file, 'r') as f, open(fname, 'w') as chord, open(fname2, 'w') as corr:
            f.readline()                                  # skip header
            chord.write('Spacer_1\tSpacer_2\tCount\n')    # write header
            corr.write('Spacer_1\tSpacer_2\tCount\n')     # write header

            for ln in f:
                cnt = ln.rstrip().split('\t')[1]
                genes = ln.rstrip().split('\t')[0].split(',')
                # we only care about pairs, skip the rest
                if len(genes) == 2:       
                    chord.write('\t'.join(genes) + '\t' + cnt + '\n')
                    corr.write('\t'.join(genes) + '\t' + cnt + '\n')
                    genes[0],genes[1] = genes[1],genes[0]
                    corr.write('\t'.join(genes) + '\t' + cnt + '\n')

        f.close()
        chord.close()
        corr.close()

def main():
   
    cmdparser = argparse.ArgumentParser(description="Count spacers + repeat in files for" 
                                       " NIH Project and produce organized files for further analysis along with different plots.",
                                        usage='%(prog)s -f <list of bam files to process> -t' 
                                        '<spacer and repeat combinations> -p <path to GitHub Repo and Rscript> [optional arguments: -d]',
                                          prog='count_spacers_NIH.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE',
                            help='File with all the bam files to process, one per line.', metavar='')
    cmdparser.add_argument('-t', '--targets',  action='store', dest='TARGETS',
                            help='spacer + repeat targets for each gene', metavar='')
    cmdparser.add_argument('-p', '--path', action='store', dest='PATH', 
                           help='Path to location of GitHub Repo and Python and R Scripts', metavar='')
    cmdparser.add_argument('-d', '--detail',  action='store_true', dest='DETAIL',
                            help='Print a more detailed description of the program.')
    cmdResults = vars(cmdparser.parse_args())
    
    cwd   = os.getcwd()     # store current working dir    
    start = time.time()     # start timer to see how long the program takes to run
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)

    if cmdResults['DETAIL']:
        print("\nmurcs_script.py")
        print("\nPurpose: Count spacers + repeat in bam files for NIH Project and ")
        print("         produce organized files for further analysis along with different plots.")
        print("\nInput: A text file with each bam file names to process and a file listing")
        print("       all the spacer+repeat sequences and names.")
        print("\n*** Please use a dedicated directory for running this pipeline.***\n")
        print("Create the directory and copy the FASTA and spacer+repeat files to the new directory.")
        print("Produce bam input file by running ls *.bam > bam_input.txt\n")
        print("Required parameters:")
        print("\t1) -f bam_input.txt")
        print("\t2) -t spacer_repeat_seq.txt")
        print("\t3) -p path_to_scripts\n")
        print("Optional arguments:")
        print("\t-d:  print a detailed description of the program.\n\n")
        print("Intermediate files are moved to the 'other_files' directory when finished")
        print("\nOriginal bam files are moved to the 'bam_files' directory when finished")
        print("\nThere are multiple outputs, each named by remove the .fasta and adding text:")
        print("\t_search_results.txt = Initial search results file")
        print("\t_gene_combinations_per_read.txt = list of genes with matches to each read")
        print("\t_gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by gene name")
        print("\t_gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries")
        print("\t  written to file (position and gene name")
        print("\t_gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted")
        print("\t  by position of match relative to the start of the read (spacer order)\n")
        print("Three output files are for all experiments combined:")
        print("\tAll_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible") 
        print("\t  constructs and sub-constructs along with counts for each experiment")
        print("\temp_sorted_possible_combinations.txt = intermediate file with all possible") 
        print("\t  construct and sub-construct combinations")
        print("\ttemp_sorted_possible_combinations_Dict.txt = intermediate file with dictionary") 
        print("\t  of all possible construct and sub-construct combinations")
        print("\t_pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with minHit of 5")
        print("\t_pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample including all hits")
        print("\t_pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample with minHit of 5")
        print("")
        print("See Kevin Myers (kmyers2@wisc.edu) with any questions.\n")
        sys.exit(1)
        
    # setup logging 
    user = pwd.getpwuid(os.getuid())[0]                # get user login name
    logging.basicConfig(filename="murcs_script-Job.log", encoding='utf-8', level=logging.INFO)
    current_time = time.ctime()
    logging.info(f' Date & time started: {current_time}')
    logging.info(f' Run by {user}')  
    
    # start timer, used to calculate total run time
    start = time.time()  

    BAM_files   = []   # hold list of initial input bam files to process
    
    if cmdResults['FILE'] is not None:
        bamfile = cmdResults['FILE']
        logging.info(f' Input bam file: "{bamfile}"')
        with open(bamfile, 'r') as f:
            for bam in f:
                BAM_files.append(bam.rstrip())
        logging.info(' List of input bam files: ')
        logging.info(BAM_files)
    else:
        print("Please provide a file with bam file names, one per line.\n")
        cmdparser.print_help()
        sys.exit(1)
    
    # create fasta files from bam file, store in list
    Fasta_files = makeFasta(BAM_files)                         
    
    logging.info(' Created the following fasta files:')
    logging.info(Fasta_files)

    # retrieve Spacer/repeat file
    if cmdResults['TARGETS'] is not None:
        gene_list = cmdResults['TARGETS']
        geneSpacerDict = geneSpacerCombinations(gene_list)
        logging.info(f' Using the following gene spacer file "{gene_list}"')
    else:
        print("Please provide a file spacer + repeat sequences.\n")
        cmdparser.print_help()
        sys.exit(1)

    # report number of fasta files to process
    number_of_files = len(Fasta_files)
    print(f"{number_of_files} Fasta files to process.\n")    
    
    # create an dictionary with all the read lengths of the original input files 
    readStats    = {}          # store a list of all lengths for all samples   
    logging.info(' Start calculating read lengths.') 
    for fsa in Fasta_files:
        seqLen, totalReads = countLines(fsa)
        if fsa not in readStats:
            readStats[fsa] = seqLen
        with open(re.sub('.fasta', '_read_lengths.txt',fsa), 'w') as out:
            for rd in seqLen:
                out.write(f'{rd}\n')
        out.close()
    logging.info(' Read lengths calculated.')
    
    # create new fasta files which only contain reads that have the repeats
    logging.info(' Filtering fasta files for reads containing repeats.')
    repeatSequence = getRepeat(gene_list)
    logging.info(' Fasta files filtered for reads which contain repeats.')
    
    logging.info(' Count repeats.')
    if repeatSequence is None:
        print('\t**** Error repeat sequence not found in spacer repeat file. ')
        logging.ERROR('**** Error repeat sequence not found in spacer repeat file. EXITING')
        cmdparser.print_help()
        sys.exit(1)
    else:
        
        argLst = []
        for fsa in Fasta_files:
            argLst.append((fsa, repeatSequence))
        argTup = tuple(argLst)
        with mp.Pool() as pool:
            res = pool.starmap(countRepeats, argTup)
    logging.info(' Count repeats complete!')
          
    # get the length of spacerRepeat
    lengthSpacerRpt = len(list(geneSpacerDict.keys())[0])

    logging.info(' Starting spacer/Repeat search')
    argLst = []
    for fasta in Fasta_files:
        argLst.append((geneSpacerDict, fasta, lengthSpacerRpt))
    argTup = tuple(argLst)
    with mp.Pool() as pool:
        res = pool.starmap(search, argTup)
    logging.info(' spacer/Repeat search complete!')
    
    #search(geneSpacerDict, fasta, lengthSpacerRpt)
    
    # count spacer/Repeats for each sample and create a combined count table
    logging.info(' Count genes')
    countDataFrame = countGenes()
    countDataFrame.to_csv('Gene_Count_Table.txt', sep="\t", index_label='spacer/Repeat', na_rep=0) 
    logging.info(' Count genes complete!')

    logging.info(' Plotting the Chord and Correlation plots.')    
    print("Plotting the Chord and Correlation plots…\n")
    # generate pairwise count tables
    makePairwiseCnt()
    chord_correlation_plots(dirPath)
    logging.info(' Chord and Correlation plotting complete!')
    
    '''
    print("Let's clean up and get out of here!\n")
    
    cleanUp( cwd )
    '''
    # end timer and do math and report how long the script took to run
    end = time.time()
    total_time = round(end - start, 2)
    total_time_min = round(total_time/60, 2)
    total_time_hours = round(total_time/60/60, 2)
    logging.info(f' Run time: {total_time_hours} hours ({total_time_min} minutes) to process the {number_of_files} FASTA files.\n')
    print(f"\nIt took {total_time_hours} hours ({total_time_min} minutes) to process the {number_of_files} FASTA files.\n")
        
if __name__ == "__main__":
    main()