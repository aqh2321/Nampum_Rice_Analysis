import os

def runcmd(name_of_file):
    os.system("samtools view -bS -1 %s.sam > %s.bam"% (name_of_file,name_of_file))
    os.system('samtools sort %s.bam -o %s.sort.bam'% (name_of_file,name_of_file))
    os.system('samtools index %s.sort.bam'% (name_of_file))
    os.system('igvtools count -z 5 -w 25 %s.sort.bam  %s.bam.tdf ref'% (name_of_file,name_of_file))
