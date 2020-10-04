
import os, sys
from time import process_time
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
time_start = process_time()
#record contain 1 referrence seq and 7 protospacer seq
#print('\n$_$finish section1: building protospacer dictionary\n')
file_name = 'Chr3 extraction_04242019.fasta'
record = list(SeqIO.parse(file_name, "fasta"))
path = sys.argv[1]
directory = os.fsencode(path)

LF = {} #LEFT FORWARD
LR = {} #LEFT REVERSE
RF = {} #RIGHT FORWARD
RR = {} #RIGHT REVERSE
origin = record[0].seq
#record contain 1 referrence seq and 7 protospacer seq

for i in range(1,8):
        spacer = record[i].seq
        spacerrev = spacer.reverse_complement()
        lenth = len(spacer)
        for m in range(0, len(origin)-lenth+1):
            comparison = origin[m:lenth]
            if spacer == comparison:
                LF[record[i].id, "LF"] = origin[(m-30):(m-20)]
                LR[record[i].id, "LR"] = origin[(m-30):(m-20)].reverse_complement()
                RF[record[i].id, "RF"] = origin[(lenth+20):(lenth+30)]
                RR[record[i].id, "RR"] = origin[(lenth+20):(lenth+30)].reverse_complement()
            elif spacerrev == comparison:
                LF[record[i].id, "LF"] = origin[(m-30):(m-20)]
                LR[record[i].id, "LR"] = origin[(m-30):(m-20)].reverse_complement()
                RF[record[i].id, "RF"] = origin[(lenth+20):(lenth+30)]
                RR[record[i].id, "RR"] = origin[(lenth+20):(lenth+30)].reverse_complement()
            else:
                m += 1
                lenth += 1
#5R and 6L are contained within a very short ~37bp seq, so replace both with 10bp in the middle, 5M
del LF['Protospacer_6', 'LF']
del RF['Protospacer_5', 'RF']
del LR['Protospacer_6', 'LR']
del RR['Protospacer_5', 'RR']
LF['Protospacer_5', 'MF'] = 'CTCGCTAACC'
LR['Protospacer_5', 'MR'] = 'GGTTAGCGAG'
print('======finish building protospacer libriary======')

'''GO TO SECTION_2'''

'''SECTION_2'''


os.system('mkdir SECTION_2')
os.system('mkdir SECTION_3')

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    sample_name = filename.rpartition('.')[0][5:]
    time_start_ = process_time()
    print(filename,sample_name)
    with open(path+filename) as read: ###have/n at the end
        with open('%s.txt'%sample_name,'w') as report:
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
    print('\n$_$ finish section2: per sequence analysis for %s' %file[:-1])
    time_finish_ = process_time()
    os.system('python -m analyze1.Section3 '+sample_name)
    print('======finish analyzing ',sample_name,'takes',time_finish_-time_start_,'======\n')
print('======start to write output files======')

os.system('python -m analyze1.Section4 '+path)
time_finish = process_time()
print('Anlysis took', time_finish-time_start)

os.system('mv *.txt SECTION_3')
