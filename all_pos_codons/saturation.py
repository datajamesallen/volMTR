"""
Protein consequence codon calculator for saturation project

Will figure out the range of possible variation from wt transcript
given only single bp mutations. This will can identify possible
DeNovo mutations.
"""

file_name = 'NM_000833.5.fasta'
cdsm = 465 # start
cdsx = 4859 # end

with open(file_name) as f:
    res = f.readlines()
    res = res[1:]
    seq = ''
    for line in res:
        seq += line.rstrip()

print(seq)

dna_to_aa = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
            'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
            'TAT':'Y','TAC':'Y','TAA':'X','TAG':'X',
            'TGT':'C','TGC':'C','TGA':'X','TGG':'W',
            'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
            'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
            'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
            'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
            'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
            'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
            'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
            'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

aa_abrev = {'Ala': 'A', 'Cys': 'C', 'Asp':'D', 'Glu':'E','Phe':'F','Gly':'G',
'His':'H','Ile':'I','Lys':'K','Leu':'L','Met':'M','Asn':'N','Pro':'P',
'Gln':'Q','Arg':'R','Ser':'S','Thr':'T','Val':'V','Trp':'W','Tyr':'Y',
'Ter':'X','del':'del','ins':'ins','dup':'dup','Del':'del','Dup':'dup','Ter':'X'}
aa_unabrev = {v: k for k, v in aa_abrev.items()}

# coding domain sequence

seq = seq.replace('\n','') # remove line char
orf = seq[cdsm-1:cdsx] # calculate the open reading frame

basepairs = ['T','C','A','G']

def codonconv(codon):
    # find the amino acid for a given codon using the dictionary dna_to_aa
    try:
        aa = dna_to_aa[codon]
        return(aa)
    except:
        return(None)

def allposs_mis_syn(codon, orig_aa):
    print(orig_aa)
    """ find all posible variants for a given codon 
    and classify as missense vs synonymous """
    mis_count = 0
    syn_count = 0
    for pos, letter in enumerate(codon):
        for mut in basepairs:
            if letter == mut:
                continue
            if pos == 0:
                mutcodon = mut + codon[1:]
            if pos == 1:
                mutcodon = codon[0] + mut + codon[2]
            if pos == 2:
                mutcodon = codon[:2] + mut
            aa = codonconv(mutcodon)
            aa = aa_unabrev[aa]
            if aa:
                if aa == 'Ter':
                    continue
                if aa == orig_aa:
                    syn_count += 1
                else:
                    mis_count += 1
            else:
                print(codon)
    return (mis_count, syn_count)

def seqconv_mis_syn(orf):
    codonlist = [['aaNum','mis_pos','syn_pos','codon','orig_aa']]
    codonorf = [orf[i:i+3] for i in range(0, len(orf), 3)]
    for i, codon in enumerate(codonorf):
        orig_aa = aa_unabrev[dna_to_aa[codon]]
        mis_count, syn_count = allposs_mis_syn(codon, orig_aa)
        codonrow = [str(i+1), str(mis_count), str(syn_count), codon, orig_aa]
        codonlist.append(codonrow)
    with open(file_name + '_codons.csv','w') as ofile:
            for row in codonlist:
                ofile.write(','.join(row)+'\n')
    print('written to file')
    return None


def allposs(codon, orig):
    # find all posible protein consequences given a single codon
    aalist = []
    for pos, letter in enumerate(codon):
        for mut in basepairs:
            if pos == 0:
                mutcodon = mut + codon[1:] 
            if pos == 1:
                mutcodon = codon[0] + mut + codon[2]
            if pos == 2:
                mutcodon = codon[:2] + mut
            aa = aa_unabrev[codonconv(mutcodon)]
            if aa == None:
                continue
            aalist.append(aa)
    aaset = set(aalist)
    if orig in aaset:
        aaset.remove(orig)
    codonlist = []
    for pos, letter in enumerate(codon):
        for mut in basepairs:
            if pos == 0:
                mutcodon = mut + codon[1:] 
            if pos == 1:
                mutcodon = codon[0] + mut + codon[2]
            if pos == 2:
                mutcodon = codon[:2] + mut
            aa = aa_unabrev[codonconv(mutcodon)]
            if aa in aaset:
                codonlist.append(mutcodon)
                aaset.remove(aa)
    codonset = set(codonlist)
    return(codonset)
def primermaker(codonstart, codonend, codonseq, i, n, mutcodon):
    # creates a primer based on the sequence
    if n == 0:
        return(codonstart + mutcodon + codonend)
    if (i-n < 0):
        codonstart = "" + codonstart
    else:
        codonstart = codonstart + (codonseq[i-n])
    try:
        codonend = (codonseq[i+n]) + codonend
    except:
        codonend = codonend + ""
    return(primermaker(codonstart, codonend, codonseq, i, n-1, mutcodon))
def reversecomp(sequence):
    # reverse compliment of the primer
    rseq = ""
    for base in sequence:
        if base == 'A':
            rseq = 'T' + rseq
        if base == 'T':
            rseq = 'A' + rseq
        if base == 'G':
            rseq = 'C' + rseq
        if base == 'C':
            rseq = 'G' + rseq
    return(rseq)
def seqconv(orf, wmode = 'p'):
    """
    figures out all possible snp mutation of the given transcript

    wmode = 'p' --> print to command line
    wmode = 'w' --> print and write to file
    """
    codonlist = []
    codonorf = [orf[i:i+3] for i in range(0, len(orf), 3)]
    for i, codon in enumerate(codonorf):
        orig_aa = aa_unabrev[dna_to_aa[codon]]
        allposs_codon = allposs(codon, orig_aa)
        for mutcodon in allposs_codon:
            mutaa = aa_unabrev[codonconv(mutcodon)]
            mut_abrev = aa_abrev[mutaa]
            orig_abrev = aa_abrev[orig_aa]
            name = orig_abrev + str(i+1) + mut_abrev
            hgvs_pro = 'p.' + orig_aa + str(i+1) + mutaa
            #primer = primermaker("","",codonorf, i, 6, mutcodon)
            #codonrow = [orig_aa, str(i+1), mutaa, codon, mutcodon, primer,
            #        reversecomp(primer)]
            codonrow = [orig_aa, str(i+1), mutaa, codon, mutcodon]
            print(codonrow)
            codonlist.append(codonrow)
    if wmode == 'w':
        with open(file_name + '_codons.csv','w') as ofile:
            for row in codonlist:
                ofile.write(','.join(row)+'\n')      
        print('written to file')
    return(None)

#seqconv(orf, 'w')
seqconv_mis_syn(orf)
