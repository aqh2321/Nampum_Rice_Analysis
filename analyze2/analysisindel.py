'''please make sure bwa and biopython are installed correctly'''
'''a demo for single file analysis is available in jupyter notebook (handle_sam.jpynb)
    This script demonstrates how indels are counted
    custumed parts:
    line 36 are indexes of protospacers;
    line 37 are specific wt sequences;
    line 129 are names of the leaves.
'''

import os,sys
import pandas as pd
import numpy as np
from functools import reduce
import re
from analyze2.tdf_processing import runcmd

#get into the folder stored all the files


path = sys.argv[1]
directory = os.fsencode(path)

perfect_match_requirement = int(sys.argv[2])

def unwanted_labels(element):
    pattern = r"(\d+)([a-zA-Z]+)"
    match = re.search(pattern, element)
    if match:
        m = re.findall(pattern,element)
        if len(m) >1 :
            #print("Full match: % s" % (match.group(0)),m)
            return element
        else:
            return 'wt'
    else:
        return False

p_bank = {'p1':238,'p2':744,'p3':941,'p4':1277,'p5':1366,'p6':1402,'p7':1497}
wt_bank = {'p1':"CACAACGACG",'p2':"GCTACGTACG",'p3':"CGGGGAATCT",'p4':"GGTGGGCCCG",'p5':'CGCGAACGCG','p6':'CCCACGTTGC','p7':'CACATGTTTT'}

acceptable_index_range = {key: {key+'start': val-perfect_match_requirement+1, key+'end': val-1}  for key, val in p_bank.items()}
specific_index_range = {key: val-5 for key, val in p_bank.items()}

def protospacer_finder(x, range_dict):
    find = False
    for key,value_ in range_dict.items():
        #print(value_.values(), type(value_.values()))
        range_to_find = list(value_.values())
        #print(range_to_find)
        if x in range(range_to_find[0], range_to_find[1]):
            find = True
            return key
        else:
            pass
    if not find:
        #target seq not at protospacer region
        return find

frames = []
for file in os.listdir(directory):
    filename = os.fsdecode(file)
    #read all seqs in a file
    #if filename.endswith("rmdup.sorted.sam"):
    if filename.endswith(".sam"):
        #leaf_name = filename[:-17]
        print(leaf_name)
        leaf_name = filename.rpartition('_')[0][4:]
        #change skiprows to 1 if use rmdup files, 2 if use original sam files
        #df = pd.read_csv('%s_rmdup.sorted.sam'%leaf_name,skiprows = 1,delim_whitespace=True,header=None, index_col=None,names =  ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','M1','M2','M3','M4','M5','M6','M7','M8','M9'])
        df = pd.read_csv('aln_%s_kmer.sam'%leaf_name,skiprows = 2,delim_whitespace=True,header=None, index_col=None,names =  ['QNAME','FLAG','RNAME','POS','MAPQ','CIGAR','RNEXT','PNEXT','TLEN','SEQ','QUAL','M1','M2','M3','M4','M5','M6','M7','M8','M9'])
        df = df.iloc[:, : 10]
        #for duplicated sam files
        #runcmd('%s_rmdup.sorted'%leaf_name)
        df['protospacer'] = 'Null'
        unique_labels = list(df.CIGAR.unique())
        possible_indel_rows = df.loc[df['CIGAR'].apply(unwanted_labels) != False]
        possible_indel_rows = possible_indel_rows.loc[possible_indel_rows['POS'].apply(protospacer_finder,args=(acceptable_index_range,)) != False]
        possible_indel_rows.reset_index(drop=True, inplace=True)
        possible_indel_rows['real_pam_seq'] = 'Null'
        for i in range(possible_indel_rows.shape[0]):
            finder_result = protospacer_finder(possible_indel_rows.iloc[i]['POS'],acceptable_index_range)
            assert_if_wt = unwanted_labels(possible_indel_rows.iloc[i]['CIGAR'])
            possible_indel_rows.at[i,'protospacer'] = finder_result
            possible_indel_rows.at[i,'CIGAR'] = assert_if_wt
            #taken the seq from the segments
            starting_PAM_index = specific_index_range[finder_result] - possible_indel_rows.iloc[i]['POS']
            if starting_PAM_index >= 0 and possible_indel_rows.iloc[i]['POS']+starting_PAM_index+10 < possible_indel_rows.iloc[i]['POS']+25:
                possible_indel_rows.at[i,'real_pam_seq'] = possible_indel_rows.iloc[i]['SEQ'][starting_PAM_index:starting_PAM_index+10]
                if assert_if_wt == 'wt':
                    if possible_indel_rows.loc[i,'real_pam_seq'] != wt_bank[possible_indel_rows.loc[i,'protospacer']]:
                        possible_indel_rows.at[i,'real_pam_seq'] = 'Null'
                    else:
                        possible_indel_rows.at[i,'real_pam_seq'] = 'wt'
        indel_count_series = possible_indel_rows.groupby(['protospacer', 'real_pam_seq']).size()
        new_df = indel_count_series.to_frame(name = 'count').reset_index()
        new_df = new_df[new_df.real_pam_seq != 'Null']
        new_df['name'] = leaf_name
        frames.append(new_df)
        print(leaf_name, 'finish creating indel dataframe')
# this is the original form of output
output = pd.concat(frames,ignore_index = True)
'''
#find the largest number in wt first, set this as a threshold to filter out cross contamination
#if you want to use this function, unmute line 99, 100 and line 112
max_count_inwt = max(output.loc[(output['name']=='WT') & (output['real_pam_seq'] != 'wt')]['count'])
output['count'].where(output['count']>=max_count_inwt,0,inplace = True)
'''
output.to_csv('indel_count_from_sam_files_ordered_by_plants.csv',index = False)
print('please wait, another format of dataframe is being generated')
# transform the sheet so that the diversity analysis pipeline can accept
reindex = output.groupby(['protospacer', 'real_pam_seq']).size()
reindex_df = reindex.to_frame(name = 'occurence').reset_index()
frames2 = []
for frame in frames:
    output_array = []
    name = frame['name'].unique()[0]
    '''
    frame['count'].where(frame['count']>max_count_inwt,0,inplace = True)
    '''
    for i in range(reindex_df.shape[0]):
        find = False
        for j in range(frame.shape[0]):
            if frame.iloc[j]['protospacer'] == reindex_df.iloc[i]['protospacer'] and frame.iloc[j]['real_pam_seq'] == reindex_df.iloc[i]['real_pam_seq'] :
                output_array.append([frame.iloc[j]['protospacer']+frame.iloc[j]['real_pam_seq'], frame.iloc[j]['count']])
                find = True
        if not find:
            output_array.append([reindex_df.iloc[i]['protospacer']+reindex_df.iloc[i]['real_pam_seq'],0])
    reformed = np.asarray(output_array)
    df_reformd = pd.DataFrame(reformed,columns=['Type','count'])
    df_reformd.rename(columns={'count': name},inplace=True)
    frames2.append(df_reformd)

df_final = reduce(lambda left,right: pd.merge(left,right,on='Type'), frames2)

titles = ['Type','WT','WT_0.25rxn','WT1A','WT2A','WT3A','1A','1A_0.25rxn','1B','1C','2A','2A_0.25rxn','2B','3A','3A_0.25rxn','3B','4B','4C','5A','5B','5C','6A','6B','7A','7B','7C','8A','8B','9A','9B','10A','10B',\
         '11A','11B','13A','13B','14A','14B','15D','16A','16B','16C','17A','17B','17C','18A','18B','19A','19B','20A','20B','20C'\
         '21A','21B']
# output for diversity analysis pipeline
df_final = df_final.reindex(columns=titles)
df_final.to_csv('indel_count_from_sam_files_ordered_by_mutation_types.csv',index = False)
os.system('mkdir samfiles')
os.system('mv *.sam samfiles')
