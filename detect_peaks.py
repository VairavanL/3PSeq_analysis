#!/usr/bin/python

import sys
from joblib import Parallel, delayed
import multiprocessing


#Detecting the peaks from the chromosome ranges

#Following is the input file for the script (Vairavan's input) 
'''
First column - Contig name
Second column - Strand 
Third column - Alignment begin
Fourth column - Alignment end
Fifth column - Read count
Sixt column - The lenght of the range #Not required

Contig22218	-	6342	6390	1	48
Contig22218	-	9626	9675	5451	49
Contig22218	-	9627	9676	55	49
Contig22218	-	9629	9677	42	48
Contig22218	-	9630	9676	13	46
Contig22218	-	9649	9698	1	49
Contig4549	-	41328	41373	20	45
Contig4549	-	41489	41536	71	47
Contig4549	-	41492	41541	205	49
Contig4549	-	41493	41542	140	49
'''
#Note: The above file need not be sorted

if len(sys.argv) < 2:
   print "Usage: ./detect_peaks <bedcount_file> <readcut-off>"
   sys.exit (1)

rangefile = open(sys.argv[1], 'r')
contig_dict = {}
peak_dict = {}

all_contigs = []

i = 0

for line in rangefile:
    i += 1
    allrecords = line.split('\t')
    contigname = allrecords[0]
    orientation = allrecords[1]
    readrangestart = allrecords[2]
    readrangeend = allrecords[3]
    count = allrecords[4]
    if int(count) > int(sys.argv[2])-1:  #### Change to get different cut-off here...by default the coverage has to me atleast 3
       all_contigs.append((contigname,orientation))
       contig_dict[i] = [str(contigname),orientation,readrangestart,readrangeend,count]

uniq_contigs = list(set(all_contigs))


#The main function to merge the co-ordinates that have mapped
#Not a novel algorithm but got the directions from the following link
#http://stackoverflow.com/questions/5679638/merging-a-list-of-time-range-tuples-that-have-overlapping-time-ranges
#It is important to add 1 (+1) to ensure that peaks differing by 1 is merged

def merge(ranges):
    saved = list(ranges[0])
    for st, en in sorted([sorted(t) for t in ranges]):
        if st <= saved[1]: ##### Remove the +1 here to prevent the merging of the reads
            saved[1] = max(saved[1], en)
        else:
            yield tuple(saved)
            saved[0] = st
            saved[1] = en
    yield tuple(saved)

def contigrange(uniq_contigs,contig_dict=contig_dict):
    contig_specific_list = []
    result = []
    for i in contig_dict:
        if contig_dict[i][0] == uniq_contigs[0] and contig_dict[i][1] == uniq_contigs[1]:
           contig_specific_list.append((int(contig_dict[i][2]),int(contig_dict[i][3])))
    final_ranges = list(merge(sorted(contig_specific_list)))
    for allrange in final_ranges:
        try:
           peakvalues = []
           for i in contig_dict:
               if contig_dict[i][0] == uniq_contigs[0] and contig_dict[i][1] == uniq_contigs[1] and int(contig_dict[i][2]) >= int(allrange[0]) and int(contig_dict[i][3]) <= int(allrange[1]):
                  peakvalues.append(int(contig_dict[i][4]))
           result.append([uniq_contigs, allrange[0], allrange[1], sum(peakvalues)])
        except:
           return "Something went wrong!"
    return result

#Parallel processing started
num_cores = int(sys.argv[2])  #### This is the argument for using the number of processors
results = Parallel(n_jobs=num_cores)(delayed(contigrange)(i) for i in uniq_contigs)
collapsed_results = sum(results, []) 

#Slicing_list http://stackoverflow.com/questions/2231663/slicing-a-list-into-a-list-of-sub-lists
group = lambda t, n: zip(*[t[i::n] for i in range(n)])
grouped_results = group(collapsed_results, 1)

for i in grouped_results:
    print i[0][0][0],i[0][0][1],i[0][1],i[0][2],i[0][3]

'''
for i in uniq_contigs:
    print contigrange(i)

for j in uniq_contigs:
    contig_specific_list = []
    for i in contig_dict:       
        if contig_dict[i][0] == j[0] and contig_dict[i][1] == j[1]:
           contig_specific_list.append((int(contig_dict[i][2]),int(contig_dict[i][3])))
    final_ranges = list(merge(sorted(contig_specific_list)))
    for allrange in final_ranges:
        try:
            peakvalues = []
            for i in contig_dict:
                if contig_dict[i][0] == j[0] and contig_dict[i][1] == j[1] and int(contig_dict[i][2]) >= int(allrange[0]) and int(contig_dict[i][3]) <= int(allrange[1]):
                   peakvalues.append(int(contig_dict[i][4]))
            print j, allrange[0], allrange[1], sum(peakvalues)
        except:
            print "Something went wrong with", j

'''
