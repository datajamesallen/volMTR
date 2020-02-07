"""
create a file with the following headers
given:
    the csv file name from gnomAD
    the PDB file to use

aaNum, AC_syn, AC_mis, x, y, z
aaNum : Amino acid number (codon number)
AC_syn : the allele count for synonymous variants at that codon
AC_miss : the allele count for missense variants at that codon
x, y, z : the coordinates for the amino acid in the PDB file

then calculate the MTR score and percentiles, and write the output file
"""

pdbpath = 'h2A-SWISS-6ira.1.A.pdb'
gnomAD_file = 'ExACv1_GRIN2A.csv'
#gnomAD_file = 'GRIN2A.csv'

import matplotlib.pyplot as plt
import numpy as np
from pandasql import sqldf
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
parser = PDBParser(PERMISSIVE=1)
import Bio
import re
import math
from datetime import datetime

startTime = datetime.now()

def get_numeric(data):
    return float(''.join(re.findall('[0-9]', data)))

gmd_df = pd.read_csv(gnomAD_file)
gmd_df['aaNum'] = gmd_df['Consequence'].apply(get_numeric)
mis = gmd_df['Annotation'] == 'missense_variant'
syn = gmd_df['Annotation'] == 'synonymous_variant'
mis_df = gmd_df[mis]
syn_df = gmd_df[syn]
mis_df = sqldf('select aaNum, sum(CASE WHEN `Allele Count` > 0 THEN 1 ELSE 0 END) as AC from mis_df group by aaNum')
syn_df = sqldf('select aaNum, sum(CASE WHEN `Allele Count` > 0 THEN 1 ELSE 0 END) as AC from syn_df group by aaNum')
agg_df = sqldf('select mis_df.aaNum, mis_df.AC as mis_AC, syn_df.AC as syn_AC from mis_df left join syn_df on syn_df.aaNum = mis_df.aaNum union select syn_df.aaNum, mis_df.AC as mis_AC, syn_df.AC as syn_AC from syn_df left join mis_df on syn_df.aaNum = mis_df.aaNum')

# incorperate expected counts

pos_df = pd.read_csv('all_pos_codons/NM_000833.5.fasta_codons.csv')
pos_df = sqldf('select aaNum, mis_pos, syn_pos from pos_df')
print(pos_df)
agg_df = sqldf('select pos_df.*, agg_df.mis_AC, agg_df.syn_AC from pos_df left join agg_df on pos_df.aaNum = agg_df.aaNum')
real_df = agg_df.copy()

print(real_df)

# replicate the real MTR score

def calc_real_MTR(row):
    """ calculates the real MTR score on a per-row basis """
    #print(row)
    temp = real_df.copy()
    window = 31
    index = int(row['aaNum']-1)
    maxindex = len(temp)
    if index <= 15:
        res = temp.iloc[0:index+16]
    elif index+16 > maxindex:
        res = temp.iloc[index-15:maxindex]
    else:
        res = temp.iloc[index-15:index+16]
    winsize = len(res)
    mis_count = res.mis_AC.sum()
    syn_count = res.syn_AC.sum()
    #print(mis_count, syn_count, sep=' ')
    mis_pos = res.mis_pos.sum()
    syn_pos = res.syn_pos.sum()
    if (mis_count + syn_count) > 0:
        MTR = ((mis_count/(mis_count + syn_count))/(mis_pos/(mis_pos + syn_pos)))
    else:
        MTR = 0
    return MTR

dfsize = len(real_df) # number of residues
mis_allele = real_df['mis_AC'].sum() # number of missense alleles in the gene
syn_allele = real_df['syn_AC'].sum() # number of synonymous alleles in the gene
total_variation = mis_allele + syn_allele # total number of variants in the gene

real_df['MTR'] = real_df.apply(calc_real_MTR, axis=1)
real_df['MTR_rank'] = real_df['MTR'].rank()
real_df['MTR_centile'] = real_df['MTR_rank'].apply(lambda x: x / dfsize)

mtr_df = pd.read_csv('ExACv1_MTR_GRIN2A.txt', sep='\t')
#mtr_df['mtr']

fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
plt.xlabel('Protein sequence position')
ax1.set_ylabel('myMTR')
ax2.set_ylabel('MTR')
ax1.plot(real_df['aaNum'],real_df['MTR'])
ax2.plot(mtr_df['Protein_position'], mtr_df['MTR'])
plt.savefig('Orig_MTR_result.png')
plt.close()

real_df.to_csv('Orig_MTR_result.csv', index=False)


# create the pdb structure data

structure = parser.get_structure('h2A', pdbpath)
model = structure[0]
chain = model['A']
pdbstruct = []
for residue in chain:
    if not is_aa(residue):
        continue
    atom = residue["CA"] # using the alpha carbon
    # uppercase protein letter maps
    resi = Bio.Data.SCOPData.protein_letters_3to1[residue.get_resname()]
    pdbrow = [residue.get_id()[1]] + atom.get_coord().tolist()
    pdbstruct.append(pdbrow)
pdb_df = pd.DataFrame(data=pdbstruct, columns = ['aaNum','x','y','z'])
all_df = pdb_df.join(agg_df.set_index('aaNum'), on='aaNum')
all_df = all_df.fillna(0)

# calculate MTR
print('calculating MTR score...')

def calc_MTR_sphere_window(row):
    """
    calculates the MTR score for each row in the data frame 
    based on the nearby residues in a fixed distance r from the 
    center of each alpha carbon
    """
    c = np.array((row['x'],row['y'],row['z'])) # our center points
    r = 10 # number of angstroms to use as our "sliding window"
    mis_count = 0
    syn_count = 0
    temp = all_df.copy()
    # np.linalg.norm will calculate the euclidean distance between the two arrays
    temp['dist'] = temp.apply(lambda x: np.linalg.norm(np.array((x['x'],x['y'],x['z']))-c), axis=1)
    temp = temp.sort_values('dist', ascending=True)
    # sort the values and then break once we have reached greater than our threshold
    for i,irow in temp.iterrows():
        if irow['dist'] <= r:
            #print(irow['dist'])
            mis_count += irow['mis_AC']
            syn_count += irow['syn_AC']
        else:
            break
    #print(syn_count, mis_count, sep=' ')
    syn_pos = row['syn_pos']
    mis_pos = row['mis_pos']
    if (mis_count + syn_count) > 0:
        MTR = ((mis_count/(mis_count + syn_count))/(mis_pos/(mis_pos + syn_pos)))
    else:
        MTR = 0
    return MTR

def calc_MTR_closest_window(row):
    """
    calculates the MTR score for each row in the data frame
    based on the closest N residues from the center of each alpha carbon
    """
    c = np.array((row['x'],row['y'],row['z'])) # our center points
    N = 31 # number of nearest neighbors for calc
    mis_count = 0
    syn_count = 0
    temp = all_df.copy()
    # np.linalg.norm will calculate the euclidean distance between the two arrays
    temp['dist'] = temp.apply(lambda x: np.linalg.norm(np.array((x['x'],x['y'],x['z']))-c), axis=1)
    temp = temp.sort_values('dist', ascending=True)
    # sort the values and then break once we have reached greater than our threshold
    count = 0
    for i,irow in temp.iterrows():
        if count>N:
            break
        mis_count += irow['mis_AC']
        syn_count += irow['syn_AC']
        count+=1
    #print(syn_count, mis_count, sep=' ')
    syn_pos = row['syn_pos']
    mis_pos = row['mis_pos']
    if (mis_count + syn_count) > 0:
        MTR = ((mis_count/(mis_count + syn_count))/(mis_pos/(mis_pos + syn_pos)))
    else:
        MTR = 0
    return MTR

from sklearn.preprocessing import MinMaxScaler

def calc_MTR_sigmoid(row):
    """
    calculates the MTR score for each row in the data frame
    with weights, where the closer residues to the center points 
    contribute more to the number of missense/synonymous variants
    """
    c = np.array((row['x'],row['y'],row['z'])) # our center points
    mis_count = 0
    syn_count = 0
    temp = all_df.copy()
    # np.linalg.norm will calculate the euclidean distance between the two arrays
    temp['dist'] = temp.apply(lambda x: np.linalg.norm(np.array((x['x'],x['y'],x['z']))-c), axis=1)
    min_d = temp['dist'].min()
    max_d = temp['dist'].max()
    temp['mtr_dist_contrib'] = temp.apply(lambda x: (1-(x['dist'] - min_d)/(max_d - min_d)), axis=1)
    #print(temp)
    #temp = temp.sort_values('mtr_dist_contrib', descending=True)
    # iterate through all values and weight the contribution to the window
    for i,irow in temp.iterrows():
        contrib = irow['mtr_dist_contrib']
        if not np.isfinite(contrib):
            contrib = 1
        mis_count += (irow['mis_AC'] * contrib)
        syn_count += (irow['syn_AC'] * contrib)
    print(syn_count, mis_count, sep=' ')
    syn_pos = row['syn_pos']
    mis_pos = row['mis_pos']
    if (mis_count + syn_count) > 0:
        MTR = ((mis_count/(mis_count + syn_count))/(mis_pos/(mis_pos + syn_pos)))
    else:
        MTR = 0
    return MTR


dfsize = len(all_df) # number of residues
mis_allele = all_df['mis_AC'].sum() # number of missense alleles in the gene
syn_allele = all_df['syn_AC'].sum() # number of synonymous alleles in the gene
total_variation = mis_allele + syn_allele # total number of variants in the gene

#expected_var_per_aa = (total_variation/dfsize)/2
#print(expected_var_per_aa)
#print(all_df)

all_df['MTR'] = all_df.apply(calc_MTR_closest_window, axis=1)
all_df['MTR_rank'] = all_df['MTR'].rank()
all_df['MTR_centile'] = all_df['MTR_rank'].apply(lambda x: x / dfsize)

mtr_df = pd.read_csv('MTR_table.csv', sep='\t')
#mtr_df['MTR']

fig, (ax1,ax2) = plt.subplots(2,1, sharex=True)
plt.xlabel('Protein sequence position')
ax1.set_ylabel('volMTR')
ax2.set_ylabel('MTR')
ax1.plot(all_df['aaNum'],all_df['MTR'])
ax2.plot(mtr_df['protein_position'], mtr_df['mtr'])
plt.savefig('fig_MTR_result.png')
plt.close()

all_df.to_csv('MTR_result.csv', index=False)

color_df = all_df[['aaNum','MTR']]
color_df.to_csv('color_df.csv', index=False, header = False)

print('DONE. time to complete:')
print(datetime.now() - startTime)
