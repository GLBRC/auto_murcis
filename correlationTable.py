#!/usr/bin/env python
from itertools import combinations_with_replacement
from itertools import combinations
import glob
import os
import re
import sys

workingDir = os.getcwd()
'''
foundTags = set()
with open('AC1_Input_NE_01_spacer_combinations_withReplacement_value_for_correlation_plots_forPlotting.txt', 'r') as f:
    f.readline()
    for line in f:
        tags = '\t'.join(line.split()[0:2])
        foundTags.add(tags)

zeroCnts = res.difference(foundTags)
print(zeroCnts)
print(len(zeroCnts))
'''

def allPairs(geneLst):
    """allPairs
    
    Using list of input gene names, generate a set of all possible pair sets.
    Return a set of all pairs, where each pair is a string, "tag1\ttag2".    
    
    Parameters
    ----------
    geneLst : str
        Text file containing the gene name and tag sequence.
        Here we only care about the gene name.
    """
    # get the initial gene tag input list
    tagLst = []
    with open(geneLst, 'r') as f:
        f.readline()         # skip header
        f.readline()         # skip repeat line
        for line in f:
            tagLst.append(line.split()[0])
            
    temp_pairs = list(combinations_with_replacement(tagLst, 2))
    all_pairs = set()                # store gene name pairs as string in a set of all pair strings
    for p in temp_pairs:           
        pairName = '\t'.join(p)      # gene1\tgene1 etc...
        all_pairs.add(pairName)
    return all_pairs       

def makePairwiseCnt(pair_set):
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

    parameters
    ----------
    pair_set : set
        List of all possible tag pairs, each pair is a tab delimited string 
        i.e. "tag1\ttag2"
    """
    for file in glob.glob( os.getcwd() + '/*-spacerRepeat-CNTs.txt'):
        fname = re.sub('-spacerRepeat-CNTs.txt', '', file) +  '_spacer_combinations_withoutReplacement_value_for_Chord_Diagrams_forPlotting.txt'
        # dictionary key = frozenset of gene pairs , this allows checking for matching without regard to order, value = count
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
            foundSet = set()
            # write genes with pair counts to file
            for genes,cnt in srCNTs.items():   
                geneNames = '\t'.join(list(genes))
                geneNamesRvsd = '\t'.join(reversed(list(genes)))
                foundSet.add(geneNames)
                foundSet.add(geneNamesRvsd)
                gene1, gene2 = ','.join(list(genes)).split(',')
                if gene1 == gene2:
                    corrOut.write(f'{gene1}\t{gene2}\t{str(cnt)}\n')                        
                else:
                    corrOut.write(f'{gene1}\t{gene2}\t{str(cnt)}\n')
                    corrOut.write(f'{gene2}\t{gene1}\t{str(cnt)}\n')
            
            missingPairs = pair_set.difference(foundSet)

            # write genes without pair counts to file
            for genePair in missingPairs:
                gene1, gene2 = genePair.split('\t')
                if gene1 == gene2:
                    corrOut.write(f'{gene1}\t{gene2}\t750000\n')
                else:
                    corrOut.write(f'{gene1}\t{gene2}\t0\n')
                    corrOut.write(f'{gene2}\t{gene1}\t0\n')
        corrOut.close()
        
def main():
    pairs = allPairs('gene_spacers_sequence-Nicole.txt')  # list of frozenSets    
    makePairwiseCnt(pairs)

if __name__ == "__main__":
    main()