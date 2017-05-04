#==============================================================================
# This code makes the peptide fasta file from the missense mutations in a MAF
#==============================================================================

# Import packages
from __future__ import division
from pandas import read_csv
from Bio import SeqIO
import re


#==============================================================================
# Load and parse proteome file
#==============================================================================

# Load
FILE_PROT = '/path/to/your/proteome_file.txt/'
prot_db = SeqIO.parse(FILE_PROT, 'fasta')

# Parse
prot_db_dict = {}
for seq_record in prot_db:
    prot_id = seq_record.id
    seq = str(seq_record.seq)
    prot_db_dict[prot_id] = seq


#==============================================================================
# Load MAF
#==============================================================================

# Load
FILE_MAF = '/path/to/your/maf.txt'
maf = read_csv(FILE_MAF, sep='\t', header=(1)) # Check the headers row


#==============================================================================
# Parse and check data in MAF
#==============================================================================

# Filter missense mutations
missense = []
for i in xrange(len(maf)):
    if maf.Variant_Classification[i] == 'Missense_Mutation':
        missense.append(i)
        

# Check mutant heterozygosity
allele_errors = []
for i in missense:
    if maf.Tumor_Seq_Allele1[i] == maf.Tumor_Seq_Allele2[i]:
        missense.remove(i)
        allele_errors.append(i)
    if maf.Reference_Allele[i] != maf.Tumor_Seq_Allele1[i]:
        missense.remove(i)
        allele_errors.append(i)

print '\nAllele errors:', len(allele_errors)


# Check hgvs errors
hgvs_errors = []
for i in missense:
    hgvs = maf.HGVSp_Short[i]
    matchObj = re.match(r'^p.\w\d+\w$', hgvs)
    if matchObj:
        pass
    else:
        missense.remove(i)
        hgvs_errors.append(i+1)

print '\nHGVS errors:', len(hgvs_errors)


# Check ENSP in protein fasta
ensp_not_in_db = []
for i in missense:
    ensp = maf.ENSP[i]
    if ensp in prot_db_dict:
        pass
    else:
        missense.remove(i)
        ensp_not_in_db.append(i+1)

print 'ENSP not in DB:', len(ensp_not_in_db)


# Check wt_aa in position in protein fasta file
wt_aa_error = []
for i in missense:
    hgvs = maf.HGVSp_Short[i]
    wt_aa = hgvs[2]
    pos = int(hgvs[3:-1])
    ensp = maf.ENSP[i]
    wt_prot = prot_db_dict[ensp]
    if wt_aa == wt_prot[pos-1]:
        pass
    else:
        missense.remove(i)
        wt_aa_error.append(i+1)
        
print 'wt_aa and protein db errors:', len(wt_aa_error)


#==============================================================================
# Create dict with WT and MT peptides
#==============================================================================

# Select lenght of peptides according to HLA class
class_i = True
class_ii = False

if class_i is True and class_ii is False:
    middle_aa = 9
elif class_i is False and class_ii is True:
    middle_aa = 15
else:
    raise ValueError('class_i or class_ii error')
    
    
# Create WT-MT dict
aminoterminal = 0
carboxi_terminal = 0
middle = 0
asterisks_peps = 0
asterisk_check = 0
missense_peps = {}
for i in missense:
    gene = maf.Hugo_Symbol[i]
    hgvs = maf.HGVSp_Short[i]
    mutation = gene+hgvs
    ensp = maf.ENSP[i]
    wt_aa = hgvs[2]
    mt_aa = hgvs[-1]
    pos = int(hgvs[3:-1])

    if pos <= (middle_aa-1): # Mutation in aminoterminal
        aminoterminal += 1
        wt_pep = prot_db_dict[ensp][:pos+(middle_aa-1)]
        if '*' in wt_pep: # Check * character in peptide
            asterisks_peps += 1
            print i, wt_pep
            pass
        else:
            mt_pep = wt_pep[:pos-1]+mt_aa+wt_pep[pos:]
            pos_pep = [j+1 for j in xrange(len(wt_pep)) if wt_pep[j] != mt_pep[j]]

            if len(pos_pep) != 1:
                raise ValueError('pos < 14 error 1', i+1)
            elif pos_pep[0] != pos:
                raise ValueError('pos < 14 error 2', i+1)
            else:
                pass

    elif pos > len(prot_db_dict[ensp])-middle_aa: # Mutation in carboxiterminal
        carboxi_terminal += 1
        wt_pep = prot_db_dict[ensp][pos-middle_aa:]
        if '*' in wt_pep: # Check * character in peptide
            asterisks_peps += 1
            print i, wt_pep
            pass
        else:
            mt_pep = wt_pep[:(middle_aa-1)]+mt_aa+wt_pep[middle_aa:]
            pos_pep = [k+1 for k in xrange(len(wt_pep)) if wt_pep[k] != mt_pep[k]]

            if len(pos_pep) != 1:
                raise ValueError('pos > len(ensp)-15 error 1', i+1)
            elif pos_pep[0] != middle_aa:
                raise ValueError('pos > len(ensp)-15 error 2', i+1)
            else:
                pass

    else: # Muattion in the middle of the protein
        middle += 1
        wt_pep = prot_db_dict[ensp][pos-middle_aa:pos+(middle_aa-1)]
        if '*' in wt_pep: # Check * character in peptide
            asterisks_peps += 1
            print i, wt_pep
            pass
        else:
            mt_pep = wt_pep[:(middle_aa-1)]+mt_aa+wt_pep[middle_aa:]
            pos_pep = [h+1 for h in xrange(len(wt_pep)) if wt_pep[h] != mt_pep[h]]

            if len(pos_pep) != 1:
                raise ValueError('pos in middle error 1', i+1)
            elif pos_pep[0] != middle_aa:
                raise ValueError('pos in middle error 2', i+1)
            else:
                pass

    if asterisk_check < asterisks_peps:
        asterisk_check += 1
    else:
        if mutation+'|mt_pos:'+str(pos_pep[0]) not in missense_peps:
            missense_peps[mutation+'|mt_pos:'+str(pos_pep[0])] = {'wt_pep': wt_pep,
                                                                  'mt_pep': mt_pep}
                                                                  

#==============================================================================
# Write fasta file
#==============================================================================

i = 0
file_name = 'OvCS_peptides%i.txt' % ((middle_aa*2)-1)
fasta_file = open(file_name, 'w')
for pep in sorted(missense_peps.keys()):
    i += 1
    fasta_file.writelines('>'+pep+'|wt|'+str(i)+'\n'+missense_peps[pep]['wt_pep']+'\n')
    i += 1
    fasta_file.writelines('>'+pep+'|mt|'+str(i)+'\n'+missense_peps[pep]['mt_pep']+'\n')
fasta_file.close()
