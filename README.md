# **Auto-MuRCiS Pipeline**

### Purpose:

This repository is an more user-friendly pipline used to analyze multiplex random CRISPR array assembly with high-throughput long-read sequence analysis (MuRCiS) by [Ellis et al](https://elifesciences.org/reviewed-preprints/86903). This pipeline is a single command and will take ccs_bam PacBio files, determine the read length, plot box plots, quantitate the number and type of unique CRISPRi arrays, and plot correlation and chord plots. 

**Note:** These scripts have been tested on Linux (CentOS7) and Mac (MacOS 12.6). Some modification may be required to run these on other operating systems. 


### Using `murcs_script.py` to analyze data
This script requires three input files:

- A text file listing the BAM files to process, on per line (`-f`).
- A two column (tab delimited) text file listing the sequences (5'-3') for the repeat and the different spacers used in the experiment (see example below and in Example directory) (`-t`).
- Complete path to these scripts and BAM files (`-p`).

Spacer/Repeat Sequence File should be a two column, tab delimited file with the first tab
being the name of the Spacer or repeat and the second column the sequence to search (`-t`):
    
    name        sequence
    repeat      GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC
    lpg0404a    GCAAGTTACGGATACTGTCTACAG
    lpg0404b    GGGTGCCTTGTTCTTCCTTATCCC


There are multiple outputs, each named by remove the .fasta and adding text:
    _search_results.txt = Initial search results file
    _gene_combinations_per_read.txt = list of genes with matches to each read
    Gene_Count_Table.txt = constructs with counts for each experiment
    _pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with 
      minHit of 5
    _pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample 
      including all hits
    _pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample 
      with minHit of 5
    
Requirements:
    
	Python 3
    Python modules:
        BioPython
        functools
        glob
        itertools
        logging
        multiprocessing
        pandas
        pyfastx  
    R
    R libraries:
        ggplot2
        reshape2
        circlize

Script must be run in the same directory as the BAM files and Spacer/Repeat Sequence File

To run the script, `cd` to the directory containing the files to analyze and run this as an example:

	`python ~/scripts/MuRCS_pipeline/murcs_script.py -f bam_files.txt -t Spacer_Repeat_sequence_File.txt -p ~/bin/MuRCS_pipeline`

### Notes

This has been tested on Linus (CentOS7) and MacOS (12.6). Other operating systems may require modifications.


