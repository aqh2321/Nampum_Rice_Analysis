'''
Section2: Per_sequence_analysis
'''
import os, sys
from Supporter import library
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

file_name = 'Chr3 extraction_04242019.fasta'
record = list(SeqIO.parse(file_name, "fasta"))

lib = library(record)
LF = lib.LF()
RF = lib.RF()
LR = lib.LR()
RR = lib.RR()

'''processing'''
file = sys.argv[1]

with open(file+'.fq') as read: ###replace with files ends in A/B/C/D
    with open(file+'.txt','w') as report:
        for title,seq,qual in FastqGeneralIterator(read):
            fqid = title
            fqlenth = len(seq)
            n = 10
            for fq in range(0, (fqlenth-9)):
                frequence = seq[fq:n]
                for key, value in LF.items():
                    if frequence == value:
                        report.write("%s ! %s%d \n" % (fqid, str(key[0][-1]+key[1]), fq))
                for key, value in LR.items():
                    if frequence == value:
                        report.write("%s ! %s%d \n" % (fqid, str(key[0][-1]+key[1]), fq))
                for key, value in RF.items():
                    if frequence == value:
                        report.write("%s ! %s%d \n" % (fqid, str(key[0][-1]+key[1]), fq))
                for key, value in RR.items():
                    if frequence == value:
                        report.write("%s ! %s%d \n" % (fqid, str(key[0][-1]+key[1]), fq))

                n += 1

print('\n$_$ finish section2: per sequence analysis\n')

'''GO TO SECTION_3 and 4'''
#os.sys('mkdir SECTION_4')
#os.sys('python Section4.py '+file)
