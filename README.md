# **MuRCiS Pipeline**

### Purpose:

This repository is associated with the coupling of multiplex random CRISPR array assembly with high-throughput long-read sequence analysis (MuRCiS). This will take ccs_bam PacBio files, determine the read length, plot box plots, quantitate the number and type of unique CRISPRi arrays, and plot correlation and chord plots.

The first part of the pipeline requires basic Linux commands, which are explained here. This will count the read lengths and remove reads lacking a match to the repeat region. The second part of the pipeline is more automated and uses `murcs_script.py` for this.

**Note:** These scripts have been tested on Linux (CentOS7) and Mac (MacOS 12.6). Some modification may be required to run these on other operating systems. 


### Determining the Read Lengths of the PacBio files

 

### Process the `repeat_match.fasta` files using `murcs_script.py`

This script will count the number of times each spacer + repeat combination appears in FASTA files

Input is a list of the spacer and repeat combinations and list of FASTA files to search

Spacer/Repeat Combination File should be a two column, tab delimited file with the first tab
being the name of the Spacer/Repeat combination and the second column the sequence to search (`-t`):
    
    ID                  Seq
    Gene_S_R_1          ATTAAGCCATGGCAGTGCAGACGATAGAGCACATAGCTAGCTATACGATAAAATCG
    Gene_S_R_2          TTCAAAACAGCATAGCTCTAAAACATTAAGCCCCAAAAAATAGATATGAGCCACAG
    
FASTA File should be a list with FASTA files to process in a single column, one per line (`-f`)

Path should be indicated as the location of the Rscript (`-p`)

There are multiple outputs, each named by remove the .fasta and adding text:
 - _search_results.txt = Initial search results file
 - _gene_combinations_per_read.txt = list of genes with matches to each read
 - _gene_combinations_per_read_SortedByGeneName.txt = list of genes sorted by gene name
 - _gene_combinations_per_read_dictionary_out.txt = dictionary of dictionaries written to file (position and gene name)
 - _gene_combinations_per_read_SortedBySpacerPosition.txt = list of genes sorted by position of match relative to the start of the read (spacer order)
 - All_Constructs_and_Subconstructs_with_Counts_All_Experiments.txt = All possible constructs and sub-constructs along with counts for each experiment
 - temp_sorted_possible_combinations.txt = intermediate file with all possible construct and sub-construct combinations
 - temp_sorted_possible_combinations_Dict.txt = intermediate file with dictionary of all possible construct and sub-construct combinations
 - _pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with minHit of 5
 - _pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample including all hits
 - _pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample with minHit of 5
    
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

Script must be run in the same directory as the FASTA files and Spacer/Repeat Combination File

To run the script, `cd` to the directory containing the files to analyze and run this as an example:

	`python ~/scripts/MuRCS_pipeline/murcs_script.py -f fasta_files.txt -t Spacer_Repeat_Combination_File.txt 
    
### Notes

This has been tested on Linus (CentOS7) and MacOS (12.6). Other operating systems may require modifications.


