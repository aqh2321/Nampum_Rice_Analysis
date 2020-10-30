'''please make sure bwa and biopython are installed correctly'''
'''This script generates sam files and tdf files'''

import os,sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from analyze2.tdf_processing import runcmd

#get into the folder stored all the files
os.system('mkdir intermediate_files')
#make genome file for later tdf file creation
'''if you want to enable tdf generation, install igvtools and find where
igv tools store genomes, and copy and replace the path here in line 16, and unmute line 74

os.system('$cut -f1,2 ref.fasta.fai > ref.chrom.sizes')
os.system('cp ref.chrom.sizes /home/nuozhang/anaconda3/share/igvtools-2.5.3-0/lib/genomes')
'''
path = sys.argv[1]
directory = os.fsencode(path)

start_site = 0
end_site = int(sys.argv[2])
avg_seq_length = int(sys.argv[3])
if not end_site:
    end_site = 25
else:
    end_site = int(end_site)
if not avg_seq_length:
    avg_seq_length = 150
else:
    avg_seq_length = int(avg_seq_length)

os.system("bwa/bwa index ref.fasta")

for file in os.listdir(directory):
    seq_bank_to_align = {}
    filename = os.fsdecode(file)
    sample_name = filename.rpartition('.')[0][5:]
    #read all seqs in a file
    print(path+filename)
    records = list(SeqIO.parse(path+filename, "fastq"))

    kmer_for_alignment = open('%s_kmer.fasta'%sample_name,'w')
    #create in silico seq bank for each sequence
    #25bp
    #data structure per leaf: {record_id: {starting index1: 25bp seq, index2: 25bp seq...}}
    for record in records:
        seq_bank_to_align[record.description] = {}
        length_of_seq = len(record)
        if length_of_seq >= avg_seq_length:
            for i in range(0, avg_seq_length, end_site):
                seq_bank_to_align[record.description][i] = record.seq[i:i+end_site]
                #seq_bank_to_align[record.description].append(record.seq[i:i+25])
        else:
            #print('length of the current sequence is shorter than expected (150bp), last %f seqs are ignored' \
            #% (150-length_of_seq))
            range_ = length_of_seq // end_site
            for i in range(0,range_ * end_site, end_site):
                seq_bank_to_align[record.description][i] = record.seq[i:i+end_site]

    #write files for alignment
    for key,value in seq_bank_to_align.items():
        new_record_id = key.replace(" ", "_")
        for key_,value_ in value.items():
            shorten_seq = SeqRecord(value_, id = new_record_id+'_'+str(key_))
            SeqIO.write(shorten_seq, kmer_for_alignment, 'fasta')
    kmer_for_alignment.close()
    print(sample_name,'finish creating kmers')

    os.system('bwa/bwa aln ref.fasta %s_kmer.fasta > %s_kmer.sai' % (sample_name,sample_name))
    os.system('bwa/bwa samse ref.fasta %s_kmer.sai %s_kmer.fasta > aln_%s_kmer.sam' % (sample_name,sample_name,sample_name))
    '''
    #generate tdfs if needed
    runcmd("aln_"+sample_name+'_kmer')
    '''

os.system('mv *.sai intermediate_files')
os.system('mv *.bam intermediate_files')
os.system('mv *.bai intermediate_files')
os.system('mv *kmer.fasta intermediate_files')
os.system('mkdir tdfs')
os.system('mv *.tdf tdfs')

# don't move sam files until analysisindel.py is run
