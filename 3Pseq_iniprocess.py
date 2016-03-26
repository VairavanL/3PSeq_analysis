#!/usr/bin/python

import sys
import string
import re
import os
import subprocess

import argparse

__author__ = 'Praveen Anand'
__date__ = '27/10/2015'

def get_args():
    '''This function parses and return arguments passed in'''
    # Wraper script description
    parser = argparse.ArgumentParser(
        description='Script is a wrapper for 3P-Seq pipeline')
    # Add arguments
    parser.add_argument(
        '-q', '--fastq', type=str, help='RAW FASTQ file of the reads', required=True)
    parser.add_argument(
        '-g', '--genome', type=str, help='Path to the complete genome file', required=False)
    parser.add_argument(
        '-c', '--config', type=str, help='Configuration file', required=False)
    parser.add_argument(
        '-o', '--output', type=str, help='Output directory to store the results', required=False)
    # Array for all arguments passed to script
    args = parser.parse_args()
    # Assign args to variables
    fastq = args.fastq
    genome = args.genome
    config = args.config
    output = args.output
    # Return all variable values
    return fastq, genome, config, output


#Getting the arguments
fastqfile, genome, config, output = get_args()

def process_3pseq(fastqfile):
    """
    This function processes the raw fastq file for 3P-Seq
    Following are the processing steps:
    Step 1 - Reverse complement
    Step 2 - Trim the terminal to retain only 2A's
    Step 3 - Exclude the sequence that contains N in them
    Step 4 - Exclude the sequence that is less than 20 nucleotide in length
    Keyword arguments:
    fastqfile - The input raw fastqfile
    """
    processed_fastq_file = open(fastqfile.replace('.fastq','_processed.fastq'), 'w')

    fastqraw = open(fastqfile)
    current_name = None

    compdna = string.maketrans("ATGC", "TACG")
    
    for i,line in enumerate(fastqraw):
        if i % 4 == 0:
           current_name = line.rstrip('\n')
        if i % 4 == 1:
           sequence = str(line.rstrip('\n'))
           revcomp = sequence.translate(compdna)[::-1]
           trimmed_seq = re.sub(r"A{2,}$", "AA", revcomp)
        if i % 4 == 2:
           second_name = str(line.rstrip('\n'))
        if i % 4 == 3:
           quality = str(line.rstrip('\n')[::-1][0:len(trimmed_seq)])
           if ((trimmed_seq[-2:] == 'AA') and (len(trimmed_seq) >= 20) and ('N' not in trimmed_seq)):
              processed_fastq_file.write(current_name+'\n')
              processed_fastq_file.write(trimmed_seq+'\n')
              processed_fastq_file.write(second_name+'\n')
              processed_fastq_file.write(quality+'\n')

process_3pseq(fastqfile)
