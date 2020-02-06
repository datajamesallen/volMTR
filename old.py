import hail as hl
ht = hl.read_table('data/gnomAD_processed.ht')

import shutil
import urllib.request as request
from contextlib import closing

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import is_aa
parser = PDBParser(PERMISSIVE=1)
import Bio

import numpy as np
import pandas as pd
import gzip
import os

def getpdb(ID, structure):
    base = 'ftp://ftp.rcsb.org/../../pub/pdb/data/structures/divided/pdb/'
    midcode = structure[1:-1]
    structcode = midcode + '/pdb' + structure + '.ent.gz'
    url = urljoin(base, structcode)
    pdbpath = 'data/protein/structures/' + ID + '_' + structure + ".pdb"
    with closing(request.urlopen(url)) as r:
        with open(pdbpath, 'wb') as f:
            shutil.copyfileobj(r, f)
    return(pdbpath)

import xmltodict

def swiss_request(ID, mode = 'json', provider = 'both', template = ''):
    """
    make requests to the swiss server givin a swissprot ID

    if mode is .json, it will return the json object.
    if mode is .pdb, it will return the path to the pdb file collected
    """
    # pre processing of url string
    if mode not in ('json','pdb'):
        raise ValueError('mode needs to be pdb or json')
    if provider not in ('both','swissmodel','pdb'):
        raise ValueError('provider needs to be one of: both, pdb, swissmodel')
    if provider == 'both':
        providerstring = ''
    else:
        providerstring = '?provider=' + provider
    if template == '':
        templatestring = ''
    else:
        templatestring = '&template=' + template
    base = 'https://swissmodel.expasy.org/repository/uniprot/'
    url = urljoin(base, ID + '.' + mode + providerstring + templatestring)
    # completed url, make request
    response = requests.get(url)
    if response:
        encoding = response.apparent_encoding
        if mode == 'pdb':
            pdb = response.content
            pdb = pdb.decode(encoding)
            pdbpath = 'data/protein/structures/' + ID + '_' + template + ".pdb"
            with open(pdbpath, 'w') as pdbfile:
                pdbfile.writelines(pdb)
            return pdbpath
        if mode == 'json':
            response = response.content
            decoded = json.loads(response)
            return decoded
    #else:
    #    return response


def whichchain(ID, template):
    """ figure out which chain belongs to the UniProt ID"""
    url = 'http://rcsb.org/pdb/rest/describeMol?structureId=' + template
    response = requests.get(url)
    if response:
        response = response.content
        o = xmltodict.parse(response)
        d = json.dumps(o)
        d = json.loads(d)
        base = d['molDescription']["structureId"]["polymer"]
        ID_matching_chains = []
        for polymer in base:
            ID_in_pdb = polymer["macroMolecule"]["accession"]["@id"]
            if ID_in_pdb == ID:
                chain = polymer["chain"]
                if isinstance(chain, list):
                    for sub_chain in chain:
                        ID_matching_chains.append(sub_chain["@id"])
                else:
                    ID_matching_chains.append(chain["@id"])
        return(ID_matching_chains)

def notlobset(ID, ht):
    fig = plt.figure()

    ht_ID = ht.filter(ht.swiss == ID)
    ht_ID = ht_ID.transmute(aa_orig = ht_ID.aa_change[0], aa_var = ht_ID.aa_change[-1])
    ht_ID = ht_ID.annotate(i = ht_ID.aa_orig + hl.str(hl.int32(ht_ID.aa_num)))
    gt = ht_ID.to_pandas()
    gt['i'] = gt['i'].astype('str')
    gt = gt.set_index('i', drop = True)
    # make a request
    result = swiss_request(ID, mode = 'json', provider = 'pdb', template = '')
    structures = []
    files = []
    chains = []
    if result:
        n_results = len(result["result"]["structures"])
        print(n_results)
        if n_results > 0:
            for n in range(0,n_results-1):
                template = result["result"]["structures"][n]["template"]
                #requests.get()
                #pdbpath = swiss_request(ID, mode = 'pdb', provider = 'pdb', template = template)
                zip_pdbpath = getpdb(ID, template)
                pdbpath = zip_pdbpath[:-2]
                pdb = open(pdbpath, "wb")
                with gzip.open(zip_pdbpath, "rb") as f:
                    bindata = f.read()
                pdb.write(bindata)
                pdb.close()
                structure = parser.get_structure(template, pdbpath)
                os.remove(pdbpath)
                match_chain = whichchain(ID, template)
                for model in structure:
                    for chain_id in match_chain:
                        ax = fig.add_subplot(111, projection='3d')
                        ax.set_aspect('equal')
                        """ this is the level at which we will run dbscan"""
                        chain = model[chain_id]
                        print(chain_id)
                        pdbstruct = []
                        for residue in chain:
                            if not is_aa(residue):
                                continue
                            atom = residue["CA"]
                            # uppercase protein letter maps
                            resi = Bio.Data.SCOPData.protein_letters_3to1[residue.get_resname()]
                            pdbrow = [resi + str(residue.get_id()[1]), residue.get_id()[1]] + atom.get_coord().tolist()
                            pdbstruct.append(pdbrow)
                        df = pd.DataFrame(data=pdbstruct, columns = ['i','aa','x','y','z'])
                        df['i'] = df['i'].astype('str')
                        df = df.set_index('i', drop = True)
                        mldf = df.join(gt) 
                        mldf = mldf.drop(['aa_num','Gene','swiss','locus.contig','locus.position','aa_orig','aa_var','alleles'], axis = 1)
                        mldf = mldf.sort_values('aa')
                        mldf['AC'] = mldf['AC'].fillna(0)
                        mldf = mldf.groupby(mldf.index).agg({'aa':np.mean,'x':np.mean,'y':np.mean,'z':np.mean,'AC':np.sum}).sort_values('aa').reset_index('i')
                        mldf = mldf.drop(['i'],axis=1)
                        mldf2 = mldf.drop(['AC'],axis=1)
                        mlmat = mldf2.as_matrix()
                        print(mlmat)
                        result = sklearn.cluster.DBSCAN().fit_predict(mlmat,sample_weight = mldf['AC'])
                        print(result)
                        
                        """
                        nogmd = mldf.loc[mldf['AC'] == 0]
                        X = nogmd['x']
                        Y = nogmd['y']
                        Z = nogmd['z']
                        gmd = mldf.loc[mldf['AC'] > 0]
                        print(gmd)
                        Xg = gmd['x']
                        Yg = gmd['y']
                        Zg = gmd['z']
                        ax.scatter(X,Y,Z,s=5,zdir='x')
                        ax.scatter(Xg,Yg,Zg,s=10,zdir='x',c='red')
                        # Create cubic bounding box to simulate equal aspect ratio
                        max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
                        Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
                        Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
                        Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
                        # Comment or uncomment following both lines to test the fake bounding box:
                        for xb, yb, zb in zip(Xb, Yb, Zb):
                           ax.plot([xb], [yb], [zb], 'w')
                        plt.grid()
                        plt.show()
                        """

                