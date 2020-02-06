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
gnomAD_file = 'GRIN2A.csv'


from pandasql import sqldf
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
parser = PDBParser(PERMISSIVE=1)
import Bio
import re
import math

def get_numeric(data):
    return float(''.join(re.findall('[0-9]', data)))

gmd_df = pd.read_csv(gnomAD_file)
gmd_df['aaNum'] = gmd_df['Consequence'].apply(get_numeric)
mis = gmd_df['Annotation'] == 'missense_variant'
syn = gmd_df['Annotation'] == 'synonymous_variant'
mis_df = gmd_df[mis]
syn_df = gmd_df[syn]
mis_df = sqldf('select aaNum, sum(`Allele Count`) as AC from mis_df group by aaNum')
syn_df = sqldf('select aaNum, sum(`Allele Count`) as AC from syn_df group by aaNum')
agg_df = sqldf('select mis_df.aaNum, mis_df.AC as mis_AC, syn_df.AC as syn_AC from mis_df left join syn_df on syn_df.aaNum = mis_df.aaNum union select syn_df.aaNum, mis_df.AC as mis_AC, syn_df.AC as syn_AC from syn_df left join mis_df on syn_df.aaNum = mis_df.aaNum')

print(agg_df)

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

def check_sphere(r,cx,cy,cz,x,y,z):
    """ checks in the points lie within the sphere with center points and r radius """
    x1 = math.pow((x-cx), 2) 
    y1 = math.pow((y-cy), 2) 
    z1 = math.pow((z-cz), 2) 
    dist = x1 + y1 + z1 # distance between the center and point x,y,z
    if dist<(r**2):
        return True # in sphere
    elif dist ==(r**2): 
        return True # on edge of sphere
    else: 
        return False # outside of sphere

def calc_MTR(row):
    """ calculates the MTR score for each row in the data frame """
    cx = row['x']
    cy = row['y']
    cz = row['z']
    r = 10 # number of angstroms to use as our "sliding window"
    mis_count = 0
    syn_count = 0
    for i,irow in all_df.iterrows():
        #print(irow['x'])
        res = check_sphere(r,cx,cy,cz,irow['x'],irow['y'],irow['z'])
        if res:
            mis_count += irow['mis_AC']
            syn_count += irow['syn_AC']
    MTR = mis_count/(mis_count + syn_count)
    print(MTR)
    return MTR

all_df['MTR'] = all_df.apply(calc_MTR, axis=1)
all_df['MTR_rank'] = all_df['MTR'].rank()
dfsize = len(all_df)
all_df['MTR_centile'] = all_df['MTR_rank'].apply(lambda x: x / dfsize)
all_df.to_csv('MTR_result.csv', index=False)

color_df = all_df[['aaNum','MTR_centile']]
color_df.to_csv('color_df.csv', index=False, header = False)
