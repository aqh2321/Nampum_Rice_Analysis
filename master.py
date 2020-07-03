'''
This script intends to analyze the rearrangement of sequence data with 7 protospacer
Before running the script, please check the helper menu with "python master.py -h" first
files are in the same folder
'''
import os,sys
from time import process_time
import argparse

parser = argparse.ArgumentParser(description='This script works for rearrangement and indel analysis in fastq files for CRISPR-caused mutations in Rice. \
Default average length for fastq files is 150bp.\
Default length for sliding window in indel analysis is 25bp. \
Specify accordingly if changes are needed.')
parser.add_argument('input_dir', type=str, help='Input dir name that contains fq files')
parser.add_argument('-b', action='store_true',help='generate bar plots if turned on')
parser.add_argument('-c', action='store_true',help='perform high coverage over fastq sequences for indel search in stage 2 analaysis')
parser.add_argument('--output_dir', type=str, help='Output dir name for counts in csv files')
parser.add_argument('--avg', default=150,type=int, help='enter the average length of fastq sequences for analysis (default is 150bp)')
parser.add_argument('--l', default=25,type=int, help='enter the length of fragments for analysis (default is 25bp)')
parser.add_argument('--m', default=1,type=int, help='which stage of analysis do you want to start (option to jump to stage 2 analaysis)')



args = parser.parse_args()
print(args.input_dir)
print(args.l)
print(args.avg)
print(args.m)
print(args.b)
path_to_raw_data = os.getcwd()+'/'+args.input_dir+'/'
if args.b:
    #add arguments to generate bar plots
    print('yes')

if args.m == 1:
    print('start stage 1 analysis')
    #os.system('mkdir fqfiles')
    os.system('python -m analyze1.Section1'+' '+ sys.argv[1])


    print('start stage 2 analysis')
    if args.c:
        os.system('python -m analyze2.autoscan2 '+path_to_raw_data+' '+str(args.l)+' '+str(args.avg))
        os.system('python -m analyze2.analysisindel2 '+os.getcwd()+' '+str(args.l))
    else:
        os.system('python -m analyze2.autoscan '+path_to_raw_data+' '+str(args.l)+' '+str(args.avg))
        os.system('python -m analyze2.analysisindel '+os.getcwd()+' '+str(args.l))

elif args.m == 2:
    print('jump to stage 2 analysis')
    if args.c:
        os.system('python -m analyze2.autoscan2 '+path_to_raw_data+' '+str(args.l)+' '+str(args.avg))
        os.system('python -m analyze2.analysisindel2 '+os.getcwd()+' '+str(args.l))
    else:
        os.system('python -m analyze2.autoscan '+path_to_raw_data+' '+str(args.l)+' '+str(args.avg))
        os.system('python -m analyze2.analysisindel '+os.getcwd()+' '+str(args.l))

#os.system('python test.py') #library, generate libs.txt
#files is doc with file names
#files = open(sys.argv[1])

#print('finish all :) check out the dataframe!')
process_to_diversity_analysis = input('start diversity analysis?[y/n]:')
if process_to_diversity_analysis == 'y':
    print('checking needed files...')
    #add Peter's diversity analysis codes
else:
    print('finish all analysis process')
