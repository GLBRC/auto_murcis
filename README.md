# **Auto-MuRCiS Pipeline**

### Purpose:

This repository is an more user-friendly pipline used to analyze multiplex random CRISPR array assembly with high-throughput long-read sequence analysis (MuRCiS) by [Ellis et al](https://elifesciences.org/reviewed-preprints/86903). This pipeline uses a Docker image and a single Bash script to analyze ccs_bam PacBio sequencing files, determine the read length, plot box plots, quantitate the number and type of unique CRISPRi arrays, and plot correlation and Chord diagram plots. 

**Note:** The user does not need to build the Docker image from the Dockerfile or the Conda environment. It is recommended that the user clone the GitHub repository in order to obtain the bash script as well as the example files, but the Docker image should be pulled from the Docker Hub. See below for specific instructions. 

### Using `run_auto_murcis.sh` to analyze data
The user must obtain the auto_murcis Docker image from Docker hub. This will ensure a consistent environment in which the program is run. To do this:

1) Download the free [Docker Desktop](https://www.docker.com/products/docker-desktop/) for their specific hardware. A free Docker account may be required. 
2) Ensure Docker Desktop is running
3) Open the terminal for command line access ([Terminal](https://support.apple.com/guide/terminal/welcome/mac) on MacOS or using [Ubuntu](https://ubuntu.com/desktop/wsl) via the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install))
4) Pull the latest auto_murcis Docker image: `docker pull auto_murcis`

After the Docker image is installed, the entire pipeline is run using a single Bash script.

This script requires three input files:

- A text file listing the BAM files to process, on per line (`-f`). An small example BAM file is provided in the Example directory.
- A two column (tab delimited) text file listing the sequences (5'-3') for the repeat and the different spacers used in the experiment (see example below and in Example directory) (`-t`).
- A two column (comma delimited) text file associating the gene IDs with a specific color to be used in the Chord diagram plots (see example below and in the Example directory) (`-c`).

The file listing the BAM files should be a single column (`-f`):

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

To run the script, navigate to the directory containing the sequence files to use as well as the three required files and type:

`bash run_auto_murcis.sh -f seq_file.txt -t target_file.txt -c color_file.txt`

This command will start the auto_murcis Docker container, copy the files needed, and process the pipeline. After the pipeline is ended, the output will be copied back to the same directory and the Docker container will be stopped and removed.

The output folder will be named `auto_murcis_output` followed by the current date. There are multiple outputs, each named by remove the .bam and adding text:

    _search_results.txt = Initial search results file
    _gene_combinations_per_read.txt = list of genes with matches to each read
    Gene_Count_Table.txt = constructs with counts for each experiment
    _pairwise_minHit5_chordDiagram.pdf = the Chord diagram for each sample with 
      minHit of 5
    _pairwise_allHits_correlationPlot.pdf = Correlation plot for each sample 
      including all hits
    _pairwise_minHit5_correlationPlot.pdf = Correlation plot for each sample 
      with minHit of 5
    
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

