#!/usr/bin/python

from itertools import product
import operator
import numpy as np
import random
import sys
import scipy.stats

from joblib import Parallel, delayed
import multiprocessing

#Ignoring the divide and invalid errors
np.seterr(divide='ignore', invalid='ignore')

#Checking the number of arguments
if len(sys.argv) < 2:
   print "Usage:./hexamer_logic.py <fastafile> <seed>"
   sys.exit(1)

#Generating all possible motifs of given length

def motifcombination(seqstring,motiflength):
    allmotifs = []
    for i in list(product(seqstring,repeat=int(motiflength))):
        allmotifs.append(''.join(i))
    return allmotifs

allmotifs =  motifcombination('ATGC',6)

#Reading in all the fasta sequences of the peaks obtained

def sequencedict(file):
    fh = open(file)
    seqrecords = {}

    for line in fh:
        if line[0] == ">":
	   header = line[1:].rstrip()
	   seq = ""
	else:
	   seqline = line.rstrip()
	   seq += seqline
        seqrecords[str(header)] = seq
    return seqrecords

allpeaks = sequencedict(sys.argv[1])

#Determining the max length of a peak
maxlengthcontig=max(allpeaks, key=lambda k: len(allpeaks[k]))
#print maxlengthcontig,allpeaks[maxlengthcontig],len(allpeaks[maxlengthcontig])

#Defining a sliding window
def slidingwindow(fseq, window_size):
    for i in xrange(len(fseq) - int(window_size) + 1):
        yield fseq[i:i+window_size]

#Generating matrix of proper size
countmatrix={}

for i in allmotifs:
    countmatrix[i]= [0] * int(len(allpeaks[maxlengthcontig])-6+1)

#Populating the matrix
'''
for key in allpeaks:
    counter = -1
    for j in slidingwindow(allpeaks[key],6):
        counter += 1
        try:
           countmatrix[j][counter] += 1
        except KeyError:
           continue

print countmatrix['AAAAAA'][0]
'''

def populate_matrix(sequence_dict,maxlengthpeak,motifs=allmotifs):
    matrix = {}
    for i in motifs:
        matrix[i] = [0] * maxlengthpeak
    
    for key in sequence_dict:
        counter = -1
        for j in slidingwindow(sequence_dict[key],6):
            counter += 1
            try:
               matrix[j][counter] += 1
            except KeyError:
               continue
    return matrix

original_matrix = populate_matrix(allpeaks,306) 


#Starting the radomization of the sequences in the dictionary
def knuth_shuffle(x,seedint):
    random.seed(int(seedint))
    stringlist = list(x)
    for i in range(len(stringlist)-1, 0, -1):
        j = random.randrange(i + 1)
        stringlist[i], stringlist[j] = stringlist[j], stringlist[i]
    shuffstring = ''.join(stringlist)
    return shuffstring

random_sequence_dict = {}


#Generating a countmatrix for random shuffling

def random_shuffle_matrix(seedint,sequence_dict=allpeaks):
    random_sequences = {}
    for i in sequence_dict:
        random_sequences[i] = knuth_shuffle(allpeaks[i],seedint)
    return populate_matrix(random_sequences,306)


#Parallelizing the generation of the random matrices
random_seeds = []
for i in range(int(sys.argv[2]),int(sys.argv[2])+100):
    random_seeds.append(i)

#Starting the parallelization
num_cores = 30
random_matrices = Parallel(n_jobs=num_cores)(delayed(random_shuffle_matrix)(seeds) for seeds in random_seeds)

def get_statistics(motif,matrix=random_matrices,matrix2=original_matrix,length=306):
    allresults = {}
    for i in range(0,length):
        random_variables = []
        for j in matrix:
            random_variables.append(j[motif][i])
        random_mean = np.mean(random_variables)
        random_sd = np.std(random_variables)
        original_val = matrix2[motif][i]
        try:
           z_score = (original_val - random_mean)/random_sd
           p_values = scipy.stats.norm.sf(abs(z_score))
        except ZeroDivisionError:
           z_score = "NaN"
           p_values = "NA"
        allresults[i] = [motif,original_val,random_mean,random_sd,z_score,p_values] 
    return allresults

'''
aaaaaa_stats = get_statistics('AAAAAA')

for i in aaaaaa_stats:
   print i+1,aaaaaa_stats[i]


for i in random_matrices:
    print i['AAAAAA'][0]
'''

##Adapted from http://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python
def bonferroni_pval(pvalues):                
    """                                                                                                   
    consistent with R - print correct_pvalues_for_multiple_testing([0.0, 0.01, 0.029, 0.03, 0.031, 0.05, 0.069, 0.07, 0.071, 0.09, 0.1]) 
    """
    n = len(pvalues)                                                                         
    new_pvalues = []
    for i in pvalues:
         if (float(i) * n) < 1.0 :
            new_pvalues.append(float(i)*n)
         else:
            new_pvalues.append('1.0') 
    return new_pvalues

all_results = Parallel(n_jobs=num_cores)(delayed(get_statistics)(motifs) for motifs in allmotifs)

all_pvalues_list = []


all_results_till_pvalues = []

for i in all_results:
    for j in i:
       if i[j][3] > 0:
          all_pvalues_list.append(i[j][5])
          all_results_till_pvalues.append([j+1,i[j]])

adjusted_pvalues_list = bonferroni_pval(all_pvalues_list)

print len(all_pvalues_list)

for i in range(0,len(all_results_till_pvalues)):
    print all_results_till_pvalues[i][0],all_results_till_pvalues[i][1][0],all_results_till_pvalues[i][1][1],all_results_till_pvalues[i][1][2],all_results_till_pvalues[i][1][3],all_results_till_pvalues[i][1][4],all_results_till_pvalues[i][1][5],adjusted_pvalues_list[i]

'''
random_AAAAAA_values = []
for i in range(1,100):
    print "randomization step " + str(i) + " completed"
    random_AAAAAA_values.append(random_shuffle_matrix(allpeaks,i))

print random_AAAAAA_values

print populate_matrix(random_sequence_dict,330)
print random_shuffle_matrix(allpeaks,200)
print random_sequence_dict[maxlengthcontig]
print allpeaks[maxlengthcontig]
'''
