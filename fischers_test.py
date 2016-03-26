#!/usr/bin/python

from itertools import product
import sys
import scipy.stats

from joblib import Parallel, delayed
import multiprocessing
import numpy as np

#Ignoring the divide and invalid errors
np.seterr(divide='ignore', invalid='ignore')

#Checking the arguments

if len(sys.argv) < 2:
   print "Usage: ./fischers_test.py <peak_fasta_file> "
   sys.exit(1)


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

def motifcombination(seqstring,motiflength):
    allmotifs = []
    for i in list(product(seqstring,repeat=int(motiflength))):
        allmotifs.append(''.join(i))
    return allmotifs

allmotifs =  motifcombination('ATGC',6)

impsignalcount = 0

for j in allpeaks:
    if 'AATAAA' in allpeaks[j]:
       impsignalcount += 1

#print impsignalcount
#print len(allpeaks)
#print len(allpeaks) - impsignalcount
'''
with_imp = 0
without_imp = 0
 
for j in allpeaks: 
    if ('AAAAAT' in allpeaks[j]) and ('AATAAA' in allpeaks[j]):
        with_imp += 1
    if 'AAAAAT' in allpeaks[j]:
        without_imp += 1

print with_imp
print impsignalcount - with_imp
print abs(without_imp - with_imp)
print (len(allpeaks) - impsignalcount) - (without_imp - with_imp)
'''

def getfischerstest(motif, dictseq = allpeaks, impsignalcount=impsignalcount):
    with_imp = 0
    without_imp = 0
    for j in dictseq:
        if (str(motif) in dictseq[j]) and ('AATAAA' in dictseq[j]):
           with_imp += 1
        if (str(motif) in dictseq[j]):
           without_imp += 1
    result = []
    result.append(motif)
    result.append(with_imp)
    result.append(impsignalcount-with_imp)
    result.append(without_imp - with_imp)
    result.append((len(allpeaks)-impsignalcount) - (without_imp - with_imp))
    oddsratio, pvalue = scipy.stats.fisher_exact([[with_imp,impsignalcount-with_imp], [without_imp-with_imp,(len(allpeaks)-impsignalcount) - (without_imp - with_imp)]])
    result.append(oddsratio)
    result.append(pvalue)
    return result

num_cores = 30

all_fischers = Parallel(n_jobs=num_cores)(delayed(getfischerstest)(motifs) for motifs in allmotifs)

#print getfischerstest('ATTTAA')

for i in all_fischers:
    print "\t".join(i)
