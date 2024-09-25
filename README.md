# **Auto-MuRCiS Pipeline**

## Table of Contents

  * [Overview of pipeline](#overview-of-pipeline)
  * [Running `run_auto_murcis.sh` to analyze data](#running-run_auto_murcissh-to-analyze-data)
  * [Required input files for the pipeline](#required-input-files-for-the-pipeline)
  * [Output of the pipeline](#output-of-the-pipeline)
  * [Docker requirements](#docker-requirements)
  * [Conda environment](#conda-environment)
  * [Contact](#contact)

### Overview of pipeline

This repository is an more user-friendly pipline used to analyze multiplex random CRISPR array assembly with high-throughput long-read sequence analysis (MuRCiS) by [Ellis et al](https://elifesciences.org/articles/86903). This pipeline uses a Docker image and a single Bash script to analyze ccs_bam PacBio sequencing files, determine the read length, plot box plots, quantitate the number and type of unique CRISPRi arrays, and plot correlation and Chord diagram plots. 

**Note:** The user does not need to build the Docker image from the Dockerfile or create the Conda environment. It is recommended that the user clone the GitHub repository in order to obtain the bash script (as well as the example files), but the Docker image should be pulled from the Docker Hub. See below for specific instructions. 

### Running `run_auto_murcis.sh` to analyze data

1) Pull the latest auto_murcis Docker image:

    `docker pull kevinmyers/auto_murcis`

2) Clone the GitHub repository using the following command. Then copy the sequencing files and other required files INTO this directory:

    `https://github.com/GLBRC/auto_murcis.git`

3) Run the script in the same directory as the sequencing files and other required files:

    `bash run_auto_murcis.sh -f seq_file.txt -t target_file.txt -c color_file.txt`

    This command will start the auto_murcis Docker container, copy the files needed, and process the pipeline. After the pipeline is ended, the output will be copied back to the same directory and the Docker container will be stopped and removed.

### Required input files for the pipeline

**Copy the sequencing files and text files needed into the cloned GitHub repository folder so the `run_auto_murcis.sh` is in the same directory with the other files**

- All ccs_bam files generated for the MuRCiS experiment
- A text file listing the BAM files to process, one per line (`-f`). A small example BAM file is provided in the Example directory.
- A two column (tab delimited) text file listing the sequences (5'-3') for the repeat and the different spacers used in the experiment (see example below and in Example directory) (`-t`).
- A two column (comma delimited) text file associating the gene IDs with a specific color to be used in the Chord diagram plots (see example below and in the Example directory) (`-c`).

The file listing the BAM files should be a single column (`-f`). Note that all the BAM files **MUST** be in the same directory with the other files and where the script will be run:

    example.bam
    example2.bam

Spacer/Repeat Sequence File should be a two column, **tab** delimited file with the first tab
being the name of the target or repeat and the second column the sequence to search (`-t`). Note that the header **IS REQUIRED**:
    
    target    spacer
    repeat    GTTTTAGAGCTATGCTGTTTTGAATGGTCCCAAAAC
    geneA     GCAAGTTACGGATACTGTCTACAG
    geneB     GGGTGCCTTGTTCTTCCTTATCCC

The Gene and Color File should be a two column, **comma** deliminted file with the first entry being the gene ID and the second entry being the [R-specific color](https://r-charts.com/colors/) to use (`-c`). It is recommended to use the [color names](https://r-charts.com/colors/):

    geneA,aquamarine
    geneB,bisque
    geneC,yellow3


### Output of the pipeline

The output folder will be named `auto_murcis_output` followed by the current date. There are multiple outputs and sub-directories to organize the files, and many will be present for each sequencing file processed.

	├── auto_murcis_output_DATE
	│   ├── summary_stats.txt = statistics about total constructs and total spacer pairs found
	│   ├── combined_read_length_for_AllReads.txt:  all read lengths for all reads for each sequencing file
	│   ├── combined_read_length_for_repeat_match.txt:  all read lengths for reads containing a match to the repeat sequence for each sequencing file
	│   ├── Gene_Count_Table.txt:  counts for all unique constructs for each sequencing file in the order in which the spacers appear in the sequencing read
	│   ├── Gene_Count_Table_sorted.txt:  counts for all unique constructs for each sequencing file sorted in an ascending manner based on gene ID for each construct
	|
	├── ├── plots
	│   ├── ├── Boxplot_combined_read_length_for_AllReads.pdf:  Box plot of read length in all reads for each sequencing file in PDF format
	│   ├── ├── Boxplot_combined_read_length_for_repeat_match.pdf:  Box plot of read length in reads containing a match to the repeat sequence for each sequencing file in PDF format
	│   ├── ├── Violin_Boxplot_combined_read_length_for_AllReads.pdf:  Violin plot of read length in all reads for each sequencing file in PDF format
	│   ├── ├── Violin_Boxplot_combined_read_length_for_repeat_match.pdf:  Violin plot of read length in reads containing a match to the repeat sequence for each sequencing file in PDF format
	│   ├── ├── *_spacer_combinations_withoutReplacement_value_for_Chord_Diagrams_forPlotting.txt:  counts for all pairwise gene combinations WITHOUT REPLACEMENT for each sequencing file for the Chord diagram plots
	│   ├── ├── *_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt:  counts for all pairwise gene combinations WITH REPLACEMENT for each sequencing file for the correlation plots
	│   ├── ├── *__pairwise_minHit5_chordDiagram.pdf:  Chord diagram plot for gene pairs with a MINIMUM COUNT OF 5 in PDF format
	│   ├── ├── *_pairwise_allHits_correlationPlot.pdf:  Correlation plot for gene pairs with a MINIMUM COUNT OF 0 in PDF format
	│   ├── ├── *_pairwise_minHit5_correlationPlot.pdf:  Correlation plot for gene pairs with a MINIMUM COUNT OF 5 in PDF format
	|
	├── ├── fasta_files
	│   ├── ├── *.fasta:  FASTA file for all reads, converted from the BAM sequencing file
	│   ├── ├── *-repeat_match.fasta:  FASTA file only for reads with a match to the repeat sequence
	|
	├── ├── readLengths
	│   ├── ├── *_read_lengths.txt:  Text file containing the lengths of all reads for a given sequence file
	│   ├── ├── *-repeats_match_read_lengths.txt:  Text file containing the lengths for only reads with a match to the repeat sequence for a given sequence file
	|
	├── ├── sampleCnts
	│   ├── ├── *-spacerRepeat-CNTs.txt:  All unique constructs found in the sequencing file and the number of times each construct is present
	|
	├── ├── other_files
	│   ├── ├── *-search_results.txt:  Intermediate file that contains the initial search results for each read in a sequencing file.
	│   ├── ├── *-gene_combinations_per_read.txt:  Intermediate file that contains the gene names in each construct found for each read in a sequencing file.
	│   ├── ├── *-gene-order_byRead.txt:  Intermediate file that contains the gene names in each construct found for each read in a sequencing file, sorted in an ascending manner based on gene ID.


### Docker requirements

1) Download the free [Docker Desktop](https://www.docker.com/products/docker-desktop/) for their specific hardware. A free Docker account may be required. 
2) Ensure Docker Desktop is running
3) Open the terminal for command line access ([Terminal](https://support.apple.com/guide/terminal/welcome/mac) on MacOS or using [Ubuntu](https://ubuntu.com/desktop/wsl) via the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install))
4) Pull the latest auto_murcis Docker image every time before running the script:  `docker pull kevinmyers/auto_murcis`

### Conda environment    

Modules built in the auto_murcis Conda environment:
    
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
        samtools
        pysam 
    R
    R libraries:
        ggplot2
        reshape2
        circlize
        tidyverse
        r-colorbrewer

### Contact

Please contact <kmyers2@wisc.edu> with any questions.

Developed by Kevin Myers and Mike Place at the University of Wisconsin-Madison.
