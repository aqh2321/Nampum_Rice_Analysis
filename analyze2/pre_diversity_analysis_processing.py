'''This script add labels to different mutation types. A demo is available in add_del_inv_label.ipynb
    The output is two sheets, one containing labeled mutation types combining rearrangement and indels
    The other contains frequencies calculated from the former data'''

import pandas as pd
import re
import numpy as np
from itertools import chain
import os,sys
path=sys.argv[1]
#Takes the excel sheet made from analysis 1, where rearrangement mutations are classified and counted.
df_rearrangement = pd.read_excel(path+'/'+'report_rearrangement.xlsx')
df_rearrangement.fillna(0,inplace = True)

def deletion_check(protospacer, pattern):
    #check if this may contain delection
    if len(pattern[0]) == 1:
        label = []
        what_to_return = False
        #this may be deletion
        #pair every two elements together for later analysis
        pattern_pairs = [''.join(x) for x in zip(pattern[:-1], pattern[1:])]
        #this is used to concatinate protospacer once the pattern is confirmed
        #pair every two protospacer number together for later analysis
        protospacer_pairs= [[x,y] for x,y in zip(protospacer[:-1], protospacer[1:])]
        for i in range(len(pattern_pairs)):
            #everytime before another round in for loop, check if the mutation type is labeled as
            #not meaningful already
            if 'N/A' in label:
                return 'None'
            if pattern_pairs[i] == 'LR' or pattern_pairs[i] == 'LM':
                #if protospacer numbers are the same
                if protospacer_pairs[i][0] == protospacer_pairs[i][1]:
                    #this will be wt label
                    label.append('p'+protospacer_pairs[i][0]+'wt')
                    if what_to_return != 'del' and what_to_return != True:
                        what_to_return = 'wt'
                elif int(protospacer_pairs[i][0]) < int(protospacer_pairs[i][1]):
                    #this will be deletion label
                    label.append('p'.join(['', *protospacer_pairs[i]])+'del')
                    what_to_return = 'del'
                else:
                    label.append('N/A')
            elif pattern_pairs[i] == 'RL':
                #if the protospacer numbers are continuous
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1])-1:
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'del' and what_to_return != True:
                        what_to_return = 'wt'
                #if protospacer numbers are the same or the numbers are not continuous
                else:
                    label.append('N/A')
            elif pattern_pairs[i] == 'MR':
                #if the protospacer numbers are continuous
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1])-1:
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'del' and what_to_return != True:
                        what_to_return = 'wt'
                elif int(protospacer_pairs[i][0])+1 < int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'del')
                    what_to_return = 'del'
                else:
                    label.append('N/A')
            else:
                #all other combinations are n/a
                label.append('N/A')
        if 'N/A' in label:
            return 'None'
        #if deletion present in the label, return only deletion label, hide wt labels
        elif what_to_return == 'del':
            del_labels= [x for x in label if 'del' in x]
            return del_labels
        else:
            return label
    #this may be inversion
    else:
        inv_result = inversion_check(protospacer, pattern)
        return inv_result

def inversion_check(protospacer_list, pattern_list):
    label = []
    what_to_return = False
    #pair every two elements together for later analysis
    pattern_pairs = [''.join(x) for x in zip(pattern_list[:-1], pattern_list[1:])]
    #pair every two protospacer number together for later analysis
    protospacer_pairs= [[x,y] for x,y in zip(protospacer_list[:-1], protospacer_list[1:])]
    for i in range(len(pattern_pairs)):
        #everytime before another round in for loop, check if the mutation type is labeled as
        #not meaningful already
        if 'N/A' in label:
            return 'None'
        if pattern_pairs[i] == 'LFLR' or pattern_pairs[i] == 'RRRF' or pattern_pairs[i] == 'RRMF' \
        or pattern_pairs[i] == 'LFMR' or pattern_pairs[i] == 'MRRF' or pattern_pairs[i] == 'MFLR':
            #if protospacer numbers are NOT the same
            if int(protospacer_pairs[i][0]) < int(protospacer_pairs[i][1]):
                label.append('p'.join(['', *protospacer_pairs[i]])+pattern_pairs[i][0]+pattern_pairs[i][2]+'inv')
                if what_to_return != 'delinv':
                    what_to_return = 'inv'
            elif int(protospacer_pairs[i][0]) > int(protospacer_pairs[i][1]):
                protospacer_pairs[i][0],protospacer_pairs[i][1] =  protospacer_pairs[i][1],protospacer_pairs[i][0]
                label.append('p'.join(['', *protospacer_pairs[i]])+pattern_pairs[i][0]+pattern_pairs[i][2]+'inv')
                if what_to_return != 'delinv':
                    what_to_return = 'inv'
            else:
                label.append('N/A')
        elif pattern_pairs[i] == 'LFRF'or pattern_pairs[i] == 'RFLF' or pattern_pairs[i]=='LFMF'or pattern_pairs[i]=='MFRF':
            #simplify the pattern
            pattern_sent_for_simplified_check = pattern_pairs[i][::2]
            if pattern_sent_for_simplified_check == 'LR' or pattern_sent_for_simplified_check == 'LM':
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'inv' and what_to_return != True and what_to_return != 'delinv' :
                        what_to_return = 'wt'
                elif int(protospacer_pairs[i][0]) < int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'del')
                    what_to_return = 'delinv'
                else:
                    label.append('N/A')
            elif pattern_sent_for_simplified_check == 'RL':
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1])-1:
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'inv' and what_to_return != True and what_to_return != 'delinv' :
                        what_to_return = 'wt'
                else:
                    label.append('N/A')
            elif pattern_sent_for_simplified_check == 'MR':
                #if the protospacer numbers are continuous
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1])-1:
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'inv' and what_to_return != True and what_to_return != 'delinv' :
                        what_to_return = 'wt'
                elif int(protospacer_pairs[i][0])+1 < int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'del')
                    what_to_return = 'delinv'
                else:
                    label.append('N/A')
        elif pattern_pairs[i] == 'LRRR' or pattern_pairs[i] == 'RRLR' or pattern_pairs[i] == 'RRMR':
            pattern_sent_for_simplified_check = pattern_pairs[i][::2]
            if pattern_sent_for_simplified_check == 'LR':
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1])+1:
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'inv' and what_to_return != True and what_to_return != 'delinv' :
                        what_to_return = 'wt'
                else:
                    label.append('N/A')
            elif pattern_sent_for_simplified_check == 'RL':
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'inv' and what_to_return != True and what_to_return != 'delinv' :
                        what_to_return = 'wt'
                elif int(protospacer_pairs[i][0]) > int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'del')
                    what_to_return = 'delinv'
                else:
                    label.append('N/A')
            elif pattern_sent_for_simplified_check == 'RM':
                if int(protospacer_pairs[i][0]) == int(protospacer_pairs[i][1])+1:
                    label.append('p'.join(['', *protospacer_pairs[i]])+'wt')
                    if what_to_return != 'inv' and what_to_return != True and what_to_return != 'delinv' :
                        what_to_return = 'wt'
                elif int(protospacer_pairs[i][0])+1 > int(protospacer_pairs[i][1]):
                    label.append('p'.join(['', *protospacer_pairs[i]])+'del')
                    what_to_return = 'delinv'
                else:
                    label.append('N/A')
        else:
            label.append('N/A')
    if 'N/A' in label:
            return 'None'
    #if inversion present in the label, return only inversion label, hide wt labels
    elif what_to_return == 'inv':
        inv_labels= [x for x in label if 'inv'in x]
        return inv_labels
    elif what_to_return == 'delinv':
        delinv_labels = [x for x in label if 'wt' not in x]
        return delinv_labels
    else:
        return label

def pattern_check(element):
    pattern_rough = r"(\d+)([a-zA-Z]+)"
    match_rough = re.findall(pattern_rough,element)
    protospacer_list = list(''.join(x[0]) for x in match_rough)
    pattern_list = list(''.join(x[1]) for x in match_rough)
    deletion_output = deletion_check(protospacer_list, pattern_list)
    return deletion_output

row_name = []
for i in range(df_rearrangement.shape[0]):
    pattern_check_output = pattern_check(df_rearrangement.loc[i][0])
    row_name.append([pattern_check_output,len(pattern_check_output)])
df_rearrangement_label = pd.DataFrame(np.array(row_name),columns = ['mut','num_of_mut_types'])


result = pd.concat([df_rearrangement, df_rearrangement_label], axis=1, sort=False)
titles = ['Type','WT','WT_0.25rxn','WT1A','WT2A','WT3A','1A','1A_0.25rxn','1B','1C','2A','2A_0.25rxn','2B','3A','3A_0.25rxn','3B','4B','4C','5A','5B','5C','6A','6B','7A','7B','7C','8A','8B','9A','9B','10A','10B',\
         '11A','11B','13A','13B','14A','14B','15D','16A','16B','16C','17A','17B','17C','18A','18B','19A','19B','20A','20B','20C'\
         '21A','21B','mut','num_of_mut_types']
df_rearrangement_with_label = result.reindex(columns=titles)
df_rearrangement_with_label[titles[1:-2]] = df_rearrangement_with_label[titles[1:-2]].astype(int)
df_rearrangement_with_label.drop(df_rearrangement_with_label[df_rearrangement_with_label.mut == 'None'].index, inplace=True)
df_rearrangement_with_label.reset_index(inplace = True,drop=True)

#if a mutation contains more than one mutation type (num_of_mut_types>0)
#repeat that row
df_repeated = df_rearrangement_with_label.loc[np.repeat(df_rearrangement_with_label.index.values,df_rearrangement_with_label.num_of_mut_types)]
df_repeated.reset_index(drop=True,inplace= True)
#if a sequence contain more than one mutation, those mutations are stored in list
#each different mutation is needed to be separated for counting
#create a separate dataframe for these mutations independently for anchoring
#df_repeated and df_separated_mut_label will have the same index for one mutation
all_mutationtypes_in_order=list(chain(*df_rearrangement_with_label['mut']))
df_separated_mut_label = pd.DataFrame(np.array(all_mutationtypes_in_order),columns = ['sep_mut'])
#combine these tow dataframes together
df_finalized_rearrangement = pd.concat([df_repeated,df_separated_mut_label],axis=1, sort=False)
df_finalized_rearrangement.drop(['mut','num_of_mut_types'],axis=1,inplace = True)

#combine the rearrangement dataframe with the indel dataframe
#only use the indel count from the not-high-coverage one
df_indel= pd.read_csv('indel_count_from_sam_files_ordered_by_mutation_types.csv')
df_indel['sep_mut'] = 'None'
for i in range(df_indel.shape[0]):
    if 'wt' not in df_indel.loc[i]['Type']:
        df_indel.at[i,'sep_mut'] = df_indel.loc[i][0][0:2]+'indel'
    else:
        df_indel.at[i,'sep_mut'] = df_indel.loc[i][0]+'true'
df_combined = pd.concat([df_finalized_rearrangement, df_indel], sort=False,ignore_index=True)
df_combined.to_csv('combined_rearrangement_indel_mutations_with_labels.csv',index= False)

#initiate a dataframe to store frequencies
df_initiated = pd.DataFrame(np.array([[0, 0, 0,0,0,0], [0, 0, 0,0,0,0], [0, 0, 0,0,0,0],[0, 0, 0,0,0,0],[0, 0, 0,0,0,0],[0, 0, 0,0,0,0],[0, 0, 0,0,0,0]]),
                   columns=['del', 'invL','invR', 'indel','unchanged','missingfreq'],index=['p1','p2','p3','p4','p5','p6','p7'],\
                           dtype =np.float)
#output = {'p1':[],'p2':[],'p3':[],'p4':[],'p5':[],'p6':[],'p7':[]}
for i in range(df_combined.shape[0]):
    if 'wttrue' in df_combined.loc[i]['sep_mut']:
        df_initiated.at[df_combined.loc[i]['sep_mut'][:2],'unchanged'] = df_combined.iloc[i][1:(df_combined.shape[1]-1)].sum()
        missing_plants =(df_combined.loc[i].values == 0).astype(int).sum() / (df_combined.shape[1]-2)
        df_initiated.at[df_combined.loc[i]['sep_mut'][:2],'missingfreq'] = missing_plants
    elif 'indel' in df_combined.loc[i]['sep_mut']:
        df_initiated.at[df_combined.loc[i]['sep_mut'][:2],'indel'] = df_initiated.loc[df_combined.loc[i]['sep_mut'][:2],'indel']+df_combined.iloc[i][1:(df_combined.shape[1]-1)].sum()
    elif 'RRinv' in df_combined.loc[i]['sep_mut'] or  'RMinv' in df_combined.loc[i]['sep_mut'] or 'MRinv' in df_combined.loc[i]['sep_mut']:
        df_initiated.at[df_combined.loc[i]['sep_mut'][:2],'invR'] = df_initiated.loc[df_combined.loc[i]['sep_mut'][:2],'invR']+df_combined.iloc[i][1:(df_combined.shape[1]-1)].sum()
    elif 'LLinv' in df_combined.loc[i]['sep_mut'] or  'LMinv' in df_combined.loc[i]['sep_mut'] or 'MLinv' in df_combined.loc[i]['sep_mut']:
        df_initiated.at[df_combined.loc[i]['sep_mut'][:2],'invL'] = df_initiated.loc[df_combined.loc[i]['sep_mut'][:2],'invL']+df_combined.iloc[i][1:(df_combined.shape[1]-1)].sum()
    elif 'del' in df_combined.loc[i]['sep_mut']:
        df_initiated.at[df_combined.loc[i]['sep_mut'][:2],'del'] = df_initiated.loc[df_combined.loc[i]['sep_mut'][:2],'del']+df_combined.iloc[i][1:(df_combined.shape[1]-1)].sum()

#calculate frequencies for diversity anslysis
'''Equation used to count total number of mutations
total mutation count+ wt count +  total* missing percentage = total<br />
total = (total mutation count+ wt count)/(1-missing percentage)'''
df_initiated['total_mutation_and_wt_count'] =df_initiated.loc[:,'del':'unchanged'].sum(axis=1)
df_initiated['total'] =df_initiated['total_mutation_and_wt_count']/(1-df_initiated['missingfreq'])
df_initiated['delfreq'] = df_initiated['del'] / df_initiated['total']
df_initiated['invLfreq'] = df_initiated['invL'] / df_initiated['total']
df_initiated['invRfreq'] = df_initiated['invR'] / df_initiated['total']
df_initiated['indelfreq'] = df_initiated['indel'] / df_initiated['total']
df_initiated['unchangedfreq'] = df_initiated['unchanged'] / df_initiated['total']

df_initiated = df_initiated.reindex(sorted(df_initiated.columns), axis=1)

df_initiated.to_csv('frequencies and counts over all plants.csv')
print('finish generating frequencies and labeled mutation sheets')
