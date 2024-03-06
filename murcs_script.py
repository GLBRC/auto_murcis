#!/usr/bin/env python
"""
murcs_script.py

This script will count the number of times each spacer + repeat combination appears in FASTA files

Notes
-----
Script makes use of Python multiprocessing to speed up processing.

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
    -search_results.txt = Initial search results file
    -gene_combinations_per_read.txt = list of genes with matches to each read
    -gene-order_byRead.txt = list of gene in order of appears in reads
    Gene_Count_Table.txt = constructs with counts for each experiment
    _pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with 
      minHit of 5
    _pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample 
      including all hits
    _pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample 
      with minHit of 5
    
Requirements
------------
    Python 3 (specifically tested with 3.10.8 )
    Python modules:
        BioPython
        Pandas
        pyfastx
        pysam
    R
    R libraries:
        ggplot2
        reshape2
        circlize

Script must be run in the same directory as the bam and Spacer/Repeat Combindation files.

Script tested on MacOS 12.6 and CentOS Linux 7(Core). 
It may require modification to run on other operating systems

@author: kmyers2@wisc.edu, with a little help from Mike Place
"""
from Bio.Seq import Seq
from Bio import SeqIO
from datetime import date
from functools import reduce
from itertools import combinations
import argparse
import glob 
import itertools                                   # for zip_longest and combinations
import logging 
import multiprocessing as mp                       # use Pool
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

def countGenes(cntFile):
    """countGenes

    Count genes (spacer/Repeat hits) for each sample.

    Returns
    -------
    pandas data frame    
    """    
    cwd = os.getcwd()
    cnt = {}   # key = spacer/Repeat value = count
    sample = os.path.basename(re.sub('-gene_combinations_per_read.txt', '', cntFile))

    with open(cntFile, 'r') as f:
        f.readline()   # skip header
        for ln in f:
            read, genes = ln.rstrip().split('\t')
            if genes not in cnt:
                cnt[genes] = 1
            else:
                cnt[genes] += 1

    df = pd.DataFrame.from_dict(cnt, orient='index', columns=[sample])
    # write individual sample counts to file
    df.to_csv( sample + '-spacerRepeat-CNTs.txt', sep='\t', index_label='spacer/Repeat', na_rep=0)
        
    return df

def ngrams(seq, n=60):
    """ngrams
    
    Create a list of ngrams using a DNA sequence string.
    N-grams in this case are a continous list of characters.

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
    -search_results.txt = Initial search results file
    -gene_combinations_per_read.txt = list of genes with matches to each read
    -gene-order_byRead.txt = list of gene in order of appears in reads
    """
    search_terms = {}
    out_searchResults = fasta.replace('.fasta', '-search_results.txt')
    out_geneCombo     = fasta.replace('.fasta', '-gene_combinations_per_read.txt')
    out_geneOrder     = fasta.replace('.fasta', '-gene-order_byRead.txt')
    
    print(f"Working on {fasta} now…\n")
    # open fasta file to process and output files 
    with open(fasta, 'r') as f, open(out_searchResults, 'w') as sr, open(out_geneCombo,'w') as cout, open(out_geneOrder, 'w') as gout:
        # write header information
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
    if not os.path.exists('other_files'):
        os.mkdir( "other_files" )
    bamDir = cwd + "/other_files/"
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("search_results.txt") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("_per_read.txt") ]
    [ os.rename( (cwd + fn), (bamDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("order_byRead.txt") ]
    # organize the original fasta files
    if not os.path.exists('fasta_files'):
        os.mkdir( "fasta_files" )
    fastaDir = cwd + "/fasta_files/"
    [ os.rename( (cwd + fn), (fastaDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".fasta") ]
    if not os.path.exists('plots'):
        os.mkdir("plots")
    plotDir = cwd + "/plots/"
    [ os.rename( (cwd + fn), (plotDir + fn) ) for fn in os.listdir(cwd) if fn.endswith(".pdf") ]
    [ os.rename( (cwd + fn), (plotDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("forPlotting.txt") ]
    if not os.path.exists('readLengths'):
        os.mkdir("readLengths")
    lenDir = cwd + "/readLengths/"
    [ os.rename( (cwd + fn), (lenDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("read_lengths.txt") ]
    if not os.path.exists('sampleCnts'):
        os.mkdir("sampleCnts")
    sampleDir = cwd + "/sampleCnts/"
    [ os.rename( (cwd + fn), (sampleDir + fn) ) for fn in os.listdir(cwd) if fn.endswith("-CNTs.txt") ]
    # remove .fai and .fxi files
    [ os.remove(cwd + fn) for fn in os.listdir(cwd) if fn.endswith(".fai")]
    [ os.remove(cwd + fn) for fn in os.listdir(cwd) if fn.endswith(".fxi")]

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
    repeatsFile = re.sub('.fasta', '-repeat_match.fasta', fsaFile)
    readLenFile = re.sub('.fasta', '-repeats_match_read_lengths.txt', fsaFile)
    sampleName  = re.sub('.fasta', '', fsaFile)
    dat = pyfastx.Fasta(fsaFile,full_index=True)    # create a fasta object
    # open file to write the reads which contain the repeats, (filter out reads w/o repeat i.e. junk)
    with open(repeatsFile, 'w') as out, open(readLenFile, 'w') as lout:
        lout.write(f'{sampleName}\n')
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
            spacerRC = str(Seq(spacer).reverse_complement())    # utilize biopython
            
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

    Generate 2 gene count tables 1) with all unidirectional pairs (without replacement) 
    and 2) each pair represented twice (with replacement), i.e.  gene1,gene2 and another
    pair gene2,gene1 with the same count for all genes in all groupings. So if we have
    lpg3000,lpg0963,lpg2223,lpg2888 we would get the following pairs:
    (lpg3000,lpg0963), (lpg3000,lpg2223),(lpg3000,lpg2888),(lpg0963,lpg2223),(lpg0963,lpg2888),
    (lpg2223,lpg2888)

    tab delimited file expected to look like: 
    spacer/Repeat   Macro_output_B2.ccs
    lpg0026b,lpg0049a,lpg0228a,lpg1917b,lpg1917a    38
    lpg0049a,lpg0026b,lpg0228a     422
    lpg0281b,lpg0228a,lpg0228b,lpg0281a     580
    lpg1658a        200

    """
    for file in glob.glob( os.getcwd() + '/*-spacerRepeat-CNTs.txt'):
        fname = re.sub('-spacerRepeat-CNTs.txt', '', file) +  '_spacer_combinations_withoutReplacement_value_for_Chord_Diagrams_forPlotting.txt'
        # dictionary key = frozenset of gene pairs , this allows checking for matching without reguard to order, value = count
        srCNTs = {}      
        with open(file, 'r') as f:
            f.readline()                          # skip header
            for line in f:
                dat = line.rstrip().split('\t')
                cnt = int(dat[1])
                genes = dat[0].split(',')
                if len(genes) > 1:                 # We only care about pairs of genes, singletons can be excluded.
                    pairs = list(combinations(genes,2))
                    frozenPairs = [frozenset(x) for x in pairs]
                    for p in frozenPairs:
                        if p not in srCNTs:
                            srCNTs[p] = cnt
                        else:
                            srCNTs[p] += cnt    
        # write to Chord diagram input file
        with open(fname, 'w') as chordOut:
            chordOut.write('Spacer_1\tSpacer_2\tCount\n')    # write header
            for genes,cnt in srCNTs.items():
                gene1, gene2 = ','.join(genes).split(',')
                chordOut.write(f'{gene1}\t{gene2}\t{str(cnt)}\n')                        
        chordOut.close()

        # write to correlation input file
        fname2 = re.sub('-spacerRepeat-CNTs.txt', '', file) + '_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt'
        with open(fname2, 'w') as corrOut:
            corrOut.write('Spacer_1\tSpacer_2\tCount\n')     # write header
            for genes,cnt in srCNTs.items():
                gene1, gene2 = ','.join(genes).split(',')
                corrOut.write(f'{gene1}\t{gene2}\t{str(cnt)}\n')                        
                corrOut.write(f'{gene2}\t{gene1}\t{str(cnt)}\n')
        corrOut.close()

def combineLengths(fileLst):
    """combineLengths

    Combine read length files for plotting into a new file.

    Parameters
    ----------
    fileLst : list
        List of files to join together  
    Returns
    -------
    df : pandas dataframe  
    """
    df = pd.concat(fileLst)
    return df

def chord_correlation_plots(path_to_R):
    """chord_correlation_plots

    Run the Rscript to plot the Chord and Correlation plots
        Chord diagrams for all samples (minHit5 only)
        Correlation charts for all samples (all and minHit5 only)
    """    
    cwd = os.getcwd()
    program = path_to_R + '/chord_correlation_plot_script.R'
    cmd = ['Rscript', program, cwd]
    subprocess.run(cmd)  

def readLengthBoxplots(path_to_R):
    """readLengthBoxplots

    Run the Rscript to plot the read length boxplots.
    
    """
    cwd = os.getcwd()
    program = path_to_R + '/plotting_boxplots_for_read_lengths.R'
    cmd = ['Rscript', program, cwd]
    subprocess.run(cmd)  

def mergeCounts(cntFile):
    """mergeCounts

    Sum up all counts for genes combination, ignoring gene order. 
    Use a frozenset as key which allows comparisons where the gene order
    doesn't matter.   
    """
    res = {}    # key = frozenset of genes value = list of counts for samples
    with open(cntFile, 'r') as f:
        sampleNames = f.readline().rstrip().split('\t')[1:]
        for line in f:
            genes = frozenset(line.split('\t')[0].split(','))
            counts = line.rstrip().split('\t')[1:]
            counts = [ float(x) for x in counts]
            if genes not in res:
                res[genes] = counts
            else:
                for i,cnt in enumerate(counts):
                    res[genes][i] += cnt
    # create a new dictionary converting the frozenset keys to strings
    cleanRes = {','.join(list(key)):value for key, value in res.items()  }
    df = pd.DataFrame.from_dict(cleanRes, orient='index', columns=sampleNames)    
    return df

def countIndividualSpacers(gene_list):
    """countIndividualSpacers
    
    Parse construct files with values and count how often each individual gene 
    spacer (gene ID) is found in any combination of spacer constructs. 
    Adapted from script by Kevin Myers.
    
    -spacerRepeat-CNTs.txt looks like:
    spacer/Repeat   AC3_Input_NE_03
    lpg0518,lpg3000 8
    lpg2552,lpg2885 256
    """
    genes = []
    with open(gene_list, 'r') as f:
        for line in f:
            genes.append(line.split()[0])
    
    genes.remove('name')       # remove unneeded items from original gene list
    genes.remove('repeat')
    
    # process all the *-spacerRepeat-CNTs.txt
    for cntFile in glob.glob('*-spacerRepeat-CNTs.txt'):
        outName = re.sub('-spacerRepeat-CNTs.txt', '_spacer_counts.txt', cntFile)
    
        gene_dict = {}        # key = gene, value = list of counts
        with open(cntFile, 'r') as f:
            for line in f:
                construct = line.split('\t')[0]          
                count = line.rstrip().split('\t')[1]
                for g in genes:                       # assign a count to a gene
                    if g in construct:
                        if g in gene_dict:
                            gene_dict[g].append(int(count))
                        else:
                            gene_dict[g]=[]
                            gene_dict[g].append(int(count))
            # sum all of a genes counts                
            out = {k: [sum(gene_dict[k])] for k in gene_dict.keys()}
            # write result
            with open(outName, 'w') as f:
                for key,val in out.items():
                    f.write(f"{key}\t{val[0]}\n")
                    
def summary():
    """summary
    
    Create a simple summary statistics file.    
    """
        
    workingDir = os.getcwd()

    with open('summary_stats.txt', 'w') as out:
        # uniq counts for each sample
        out.write('Unique Counts\n')
        out.write('sample\t>=1 hits\t>=5 hits\ttotal\n')
        for infile in sorted(list(glob.glob('*-spacerRepeat-CNTs.txt'))):
            counts = {'oneHit':0, 'fiveHit': 0, 'total':0}
            with open(infile, 'r') as uniq:
                uniq.readline()        # skip header
                for line in uniq:
                    counts['total'] += 1
                    hits = int(line.rstrip().split()[-1])
                    if hits >= 1 and hits >= 5:
                        counts['oneHit'] += 1
                        counts['fiveHit'] += 1
                    elif hits >= 1 and hits < 5:
                        counts['oneHit'] += 1
            
            outLine = f"{re.sub('-spacerRepeat-CNTs.txt', '', infile)}\t{counts['oneHit']}\t{counts['fiveHit']}\t{counts['total']}\n"
            out.write(outLine)

        out.write("\n\n")
        
        # pairwise construct numbers
        pattern = '_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt'
        out.write('Pairwise Construct Counts\n')
        out.write('sample\t>=1 hits\t\t>=5 hits\ttotal\n')    
        for infile in sorted(list(glob.glob('*_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt'))):
            counts = {'oneHit': 0, 'fiveHit': 0, 'total': 0}
            with open(infile, 'r') as pairs:
                pairs.readline()    # skip header
                for line in pairs:
                    counts['total'] += 1
                    hits = int(line.rstrip().split()[-1])
                    
                    if hits >= 1 and hits >= 5:
                        counts['oneHit'] += 1
                        counts['fiveHit'] += 1
                    elif hits >= 1 and hits < 5:
                        counts['oneHit'] += 1
            
            outLine = f"{re.sub(pattern, '', infile)}\t{counts['oneHit']}\t{counts['fiveHit']}\t{counts['total']}\n"
            out.write(outLine)

        out.write("\n\n")
        
        # reads with repeats, total reads
        out.write('Reads with Repeats\n')
        out.write('sample\tReads w/ Repeat\ttotal\n')     
        counts = {}
        
        for infile in glob.glob('*-repeats_match_read_lengths.txt'):
            name = re.sub('-repeats_match_read_lengths.txt' ,'' ,infile)
            otherFile = name + '_read_lengths.txt'
            
            if not name in counts:
                counts[name] = {'matched':0, 'total':0}
                    
            with open(otherFile, 'r') as f:
                    counts[name]['total'] = str(len(f.readlines()) - 1)
            
            with open(infile, 'r') as f:
                    counts[name]['matched'] = str(len(f.readlines()) - 1)
        
        sortNames = sorted(list(counts.keys()))
        for n in sortNames:
            outLine = f"{n}\t{counts[n]['matched']}\t{counts[n]['total']}\n"
            out.write(outLine)            

def main():
    cmdparser = argparse.ArgumentParser(description="Count spacers + repeat in files for" 
                                       " NIH Project and produce organized files for further analysis along with different plots.",
                                        usage='%(prog)s -f <list of bam files to process> -t' 
                                        '<spacer and repeat combinations> [optional arguments: -d]',
                                          prog='count_spacers_NIH.py'  )                                  
    cmdparser.add_argument('-f', '--file',  action='store', dest='FILE',
                            help='File with all the bam files to process, one per line.', metavar='')
    cmdparser.add_argument('-t', '--targets',  action='store', dest='TARGETS',
                            help='spacer + repeat targets for each gene', metavar='')
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
        print("Optional arguments:")
        print("\t-d:  print a detailed description of the program.\n\n")
        print("Intermediate files are moved to the 'other_files' directory when finished")
        print("\nOriginal bam files are moved to the 'bam_files' directory when finished")
        print("\nThere are multiple outputs, each named by removing the .fasta and adding text:")
        print("\t-search_results.txt = Initial search results file")
        print("\t-gene_combinations_per_read.txt = list of genes with matches to each read")
        print("\t-gene-order_byRead.txt = list of gene in order of appears in reads")
        print("\tGene_Count_Table.txt = constructs with counts for each experiment")
        print("Three output files are for all experiments combined:")
        print("\tAll_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible") 
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
            header = re.sub('.fasta','',fsa)
            out.write(f'{header}\n')
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
    
    # countGenes in each sample 
    logging.info(' Count genes')
    argLst = []
    for file in glob.glob(cwd + '/*-gene_combinations_per_read.txt'):
        argLst.append([file])
    argTup = tuple(argLst)
    print(argTup)
    
    dfLst = []                 # list of data frames to write as individual tables and merge into a single table
    with mp.Pool() as pool:
        for result in pool.starmap(countGenes, argTup):
            dfLst.append(result)

    # count spacer/Repeats for each sample and create a combined count table, by individual spacer/Repeats
    countDataFrame = reduce(lambda left,right: pd.merge(left,right, how="outer", left_index=True, right_index=True),dfLst)
    countDataFrame.sort_index(axis=0, inplace=True, )
    countDataFrame.to_csv('Gene_Count_Table.txt', sep="\t", index_label='spacer/Repeat', na_rep=0)     
    logging.info(' Count genes complete!')

    # Write merged gene counts, where the order of the genes is not important, just gene presence.
    finalCnt = mergeCounts('Gene_Count_Table.txt')
    finalCnt.to_csv('Gene_Count_Table_Merged.txt', sep="\t", index_label='spacer/Repeat', na_rep=0)
    
    # Sort the Gene_Count_Table.txt genes names for user, using sorted (not ideal)
    with open('Gene_Count_Table_Merged.txt', 'r') as f, open('Gene_Count_Table_sorted.txt', 'w') as out:
            out.write(f.readline())             # write header
            for line in f:
                d = line.rstrip().split('\t')
                d[0] = ','.join(sorted([ n for n in d[0].split(',')]))               
                out.write('\t'.join(d) + '\n')         # write sorted gene names with values
    f.close()
    out.close()
    
    os.remove('Gene_Count_Table_Merged.txt')
    
    logging.info(' Plotting the Chord and Correlation plots.')    
    print("Plotting the Chord and Correlation plots…\n")
    # generate pairwise count tables
    makePairwiseCnt()
    chord_correlation_plots(dirPath)
    logging.info(' Chord and Correlation plotting complete!')
    
    # start collecting the read length files so we can create tables for plotting
    readLst = []       # all lengths included
    matchReadLst = []  # only matched read lengths included
    for lenFile in glob.glob('*_read_lengths.txt'): 
        if lenFile.endswith('repeats_match_read_lengths.txt'):
            matchReadLst.append(lenFile)
        else:
            readLst.append(lenFile)

    # create the files for plotting lengths
    # all lengths first
    plotDFs = []
    for f in readLst:
        plotDFs.append(pd.read_csv(f))
    plotTable = pd.concat(plotDFs, axis=1)
    plotTable.to_csv('combined_read_length_for_AllReads.txt', sep='\t', index=False)

    # Now matched reads only
    plotDFs = []
    for f in matchReadLst:
        plotDFs.append(pd.read_csv(f))
    plotTable = pd.concat(plotDFs, axis=1)
    plotTable.to_csv('combined_read_length_for_repeat_match.txt', sep='\t', index=False)    

    logging.info(' Generating read length box plots.')
    readLengthBoxplots(dirPath)
    logging.info(' Read length box plots complete!')
    
    logging.info(' Start counting individual spacers.')
    countIndividualSpacers(gene_list)
    logging.info(' Count individual spacers complete!')
    
    logging.info(' Create summary stats.')
    summary()
    
    logging.info(" Running clean up step.\n")
    cleanUp( cwd )
    
    # end timer and do math and report how long the script took to run
    end = time.time()
    total_time = round(end - start, 2)
    total_time_min = round(total_time/60, 2)
    total_time_hours = round(total_time/60/60, 2)
    logging.info(f' Run time: {total_time_hours} hours ({total_time_min} minutes) to process the {number_of_files} FASTA files.\n')
    print(f"\nIt took {total_time_hours} hours ({total_time_min} minutes) to process the {number_of_files} FASTA files.\n")
    
if __name__ == "__main__":
    main()