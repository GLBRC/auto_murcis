#!/usr/bin/env Rscript

#Making Boxplots for read lengths for NIH CRISPR Project

#generate read lengths from bam files in terminal:
#/opt/bifxapps/samtools-1.9/bin/samtools view NE_01.ccs.bam | awk '{print length($10)}' > NE_01_read_lengths.txt
#/opt/bifxapps/samtools-1.9/bin/samtools view NE_02.ccs.bam | awk '{print length($10)}' > NE_02_read_lengths.txt
#/opt/bifxapps/samtools-1.9/bin/samtools view NE_03.ccs.bam | awk '{print length($10)}' > NE_03_read_lengths.txt
#/opt/bifxapps/samtools-1.9/bin/samtools view NE_04.ccs.bam | awk '{print length($10)}' > NE_04_read_lengths.txt
#/opt/bifxapps/samtools-1.9/bin/samtools view NE_05.ccs.bam | awk '{print length($10)}' > NE_05_read_lengths.txt
#/opt/bifxapps/samtools-1.9/bin/samtools view NE_06.ccs.bam | awk '{print length($10)}' > NE_06_read_lengths.txt

#Import files and draw boxplot

library(tidyverse)
library(RColorBrewer)
library(ggplot2)

args <- commandArgs(TRUE)
working.dir <- args[1]
setwd(working.dir)

allRead_file <- read.table(file = 'combined_read_length_for_AllReads.txt', sep = "\t", header = TRUE)  
colNum <- length(colnames(allRead_file))
palette(brewer.pal(n = colNum, name = "Dark2"))
allRead_long_df <- allRead_file %>% gather(Key, Value)

pdf("Violin_Boxplot_combined_read_length_for_AllReads.pdf", width = 12, height = 12)

Violin_Boxplot_allReads <- 
  allRead_long_df %>%
  ggplot(aes(x=Key, y=Value, fill=Key)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  labs(x="Sample", y="Read Length (bp)") +
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(size=12), axis.text.x = element_text(angle = 90)
        ) +
  ggtitle("Violin and Boxplot of All Read Lengths") +
  scale_fill_brewer(palette="Dark2")
plot(Violin_Boxplot_allReads)
dev.off()

pdf("Boxplot_combined_read_length_for_AllReads.pdf", width = 12, height = 12)

Boxplot_allReads <- 
  allRead_long_df %>%
  ggplot(aes(x=Key, y=Value, fill=Key)) +
  geom_boxplot() +
  labs(x="Sample", y="Read Length (bp)") +
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(size=12), axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("Boxplot of All Read Lengths") +
  scale_fill_brewer(palette="Dark2")
plot(Boxplot_allReads)
dev.off()

RepeatMatchRead_file <- read.table(file = 'combined_read_length_for_repeat_match.txt', sep = "\t", header = TRUE)  
colNum <- length(colnames(RepeatMatchRead_file))
palette(brewer.pal(n = colNum, name = "Dark2"))
RepeatMatchRead_long_df <- RepeatMatchRead_file %>% gather(Key, Value)

pdf("Violin_Boxplot_combined_read_length_for_repeat_match.pdf", width = 12, height = 12)

Violin_Boxplot_repeatedMatchReads <- 
  RepeatMatchRead_long_df %>%
  ggplot(aes(x=Key, y=Value, fill=Key)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="black", alpha=0.2) +
  labs(x="Sample", y="Read Length (bp)") +
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(size=12), axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("Violin and Boxplot of Repeat Matched Read Lengths") +
  scale_fill_brewer(palette="Dark2")
print(Violin_Boxplot_repeatedMatchReads)
dev.off()

pdf("Boxplot_combined_read_length_for_repeat_match.pdf", width = 12, height = 12)

Violin_Boxplot_repeatedMatchReads <- 
  RepeatMatchRead_long_df %>%
  ggplot(aes(x=Key, y=Value, fill=Key)) +
  geom_boxplot() +
  labs(x="Sample", y="Read Length (bp)") +
  theme_classic() +
  theme(legend.position="none", plot.title = element_text(size=12), axis.text.x = element_text(angle = 90)
  ) +
  ggtitle("Boxplot of Repeat Matched Read Lengths") +
  scale_fill_brewer(palette="Dark2")
print(Violin_Boxplot_repeatedMatchReads)
dev.off()