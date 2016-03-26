#!/usr/bin/python

import sys
import string
import re
import os
import subprocess

from ConfigParser import SafeConfigParser

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
        '-q', '--fastq', type=str, help='Processed FASTQ file of the reads [output from 3Pseq_iniprocess.py]', required=False)
    parser.add_argument(
        '-g', '--genome', type=str, help='Path to the bowtie genome index file', required=True)
    parser.add_argument(
        '-c', '--config', type=str, help='Configuration file', required=True)
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
processedfastqfile, genome, config, output = get_args()

parser = SafeConfigParser()
parser.read(config)


if os.path.exists(parser.get('bowtie_options','bowtiepath')):
   print "Found bowtie installation!!"
else :
   print "Sorry bowtie not found on your system. Please enter the correct path in the config file"
   sys.exit(1)

if os.path.exists(parser.get('samtools_options','sampath')):
   print "Found samtools installation!!"
else :
   print "Sorry samtools installation not found on your system. Please enter the correct path in the config file"
   sys.exit(1)

if os.path.exists(parser.get('bedtools_options','bambedpath')):
   print "Found bam2bed!!!"
else :
   print "Sorry could not fine 'bamToBed' program...please specify a proper path.."
   sys.exit(1)


def run_bowtie3P(processedfastqfile, configfile=config):
    bowtie_cmd = "("+parser.get('bowtie_options','bowtiepath')+" -q -m "+parser.get('bowtie_options','m')+" -v " +parser.get('bowtie_options','v') +" -p "+ parser.get('bowtie_options','p')+" --chunkmbs 200 --un "+ processedfastqfile.replace('.fastq','_unmapped.fastq') +" "+ str(genome) +" "+processedfastqfile + " --sam " + processedfastqfile.replace('.fastq','_aligned.sam')+")"
    p = subprocess.Popen(bowtie_cmd,stdin=None,stdout=None, shell=True)
    p.wait()
    if p.returncode != 0:
       return "Sorry!!..Something went wrong with the bowtie run...Ensure that all the variables in the bowtie are specified properly in config file.."

    processed_sam_file = open(processedfastqfile.replace('.fastq','_filtered.sam'), 'w')

    with open(processedfastqfile.replace('.fastq','_aligned.sam')) as infile:
         for line in infile:
             if line[0] == "@":
                processed_sam_file.write(line)
             if line.split('\t')[1] == "0":
                MDZtag = str(line.split('\t')[12]).rstrip('\n')                
                if re.search(r"[A-Z]0$", MDZtag):
                   processed_sam_file.write(line)
             if line.split('\t')[1] == "16":
                MDZtag = str(line.split('\t')[12]).rstrip('\n')                
                if re.search(r"MD:Z:0[A-Z]", MDZtag):                   
                   processed_sam_file.write(line)
    return "Done with bowtie alignment and processing of sam file"

def gettingbed_file(filtered_samfile, configfile=config):
    samtools_cmd = "("+parser.get('samtools_options','sampath')+" view -bS "+ filtered_samfile + " > "+ filtered_samfile.replace('.sam','.bam')+")"
    samtoolrun = subprocess.Popen(samtools_cmd,stdin=None,stdout=None, shell=True)
    samtoolrun.wait()
    bedtools_cmd = "("+parser.get('bedtools_options','bambedpath')+" bamtobed -i "+ filtered_samfile.replace('.sam','.bam') + " > "+ filtered_samfile.replace('.sam','.bed')+ ")"
    bedtools_run = subprocess.Popen(bedtools_cmd,stdin=None,stdout=None, shell=True)
    bedtools_run.wait()
    bedfile = open(filtered_samfile.replace('.sam','.bed'))
    range_count = {}
    for line in bedfile:
        bedcolumns = line.split('\t')
        uniqrange = bedcolumns[0] +"_"+ bedcolumns[5].rstrip('\n') + "_" + bedcolumns[1] + "_" + bedcolumns[2]
        try :
          range_count[uniqrange] += 1
        except KeyError:
          range_count[uniqrange] = 1

    processed_bed_file = open(filtered_samfile.replace('.sam','.bedcount'), 'w')

    for i in range_count:
        new_list = i.split("_")
        new_list.append(str(range_count[i]))
        processed_bed_file.write("\t".join(new_list)+"\n")
    
print run_bowtie3P(processedfastqfile)
print gettingbed_file(processedfastqfile.replace('.fastq','_filtered.sam'))
