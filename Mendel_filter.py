#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on:	Fri 11/09/2020
@author: 	Joanna J. Ilska

The following script performs a Mendelian check on the genotypes provided in plink .raw format. 
For full description, see README.txt 
"""

import pandas as pd
import scipy.stats as stats
pd.options.mode.chained_assignment = None
import argparse

def Mendel(row, p, n_snps, parent1, parent2):
    '''
    Based on the parental genotypes, calculate the expected and observed genotype frequencies and decide whether to accept or reject the null hypothesis.
    '''
    
    ## Specify the cross type
    p1=row[parent1]
    p2=row[parent2]

    gens=[p1,p2]
    gens.sort()

    ## Both parents homozygous for the same allele:
    if (gens==[0.0,0.0]) | (gens==[2.0,2.0]):
        cross='zzxzz'
    ## Both parents homozygous but for the opposite alleles
    elif gens==[0.0, 2.0]:
        cross='aaxbb'
    ## One parent homo, the other hetero
    elif (gens==[0.0,1.0]) | (gens==[1.0,2.0]):
        cross='nnxnp'
    ## Both parents heterozygous
    elif gens==[1.0,1.0]:
        cross='hkxhk'

    ## Dictionary of expected frequencies
    exp={
        'hkxhk':{'hh':0.25,'hk':0.5,'kk':0.25},
        'nnxnp':{'hh':0.5,'hk':0.5},
        'zzxzz':{'hh':1},
        'aaxbb':{'hk':1}}

    ## Extract the genotypes
    genotypes=row.iloc[2:]

    ## Remove missing values
    test_row=genotypes.dropna().reset_index()
    test_row.columns=['ID','genotype']

    ## Replace the genotype labels - I will just use hh, hk and kk format here. 
    test_row.replace(to_replace=0.0, value='hh', inplace=True)
    test_row.replace(to_replace=1.0, value='hk', inplace=True)
    test_row.replace(to_replace=2.0, value='kk', inplace=True)


    ## Count non-missing genotypes
    nonmis=test_row.shape[0]

    ## Count observed genotypes and put into a dictionary
    obs=test_row['genotype'].value_counts().to_dict()


    ## Remove observed genotypes not possible for a given cross.
    obs_keys_invalid=set(obs.keys())-set(exp[cross].keys())

    len(obs_keys_invalid)
    if len(obs_keys_invalid)>0:
        for x in obs_keys_invalid:
            del obs[x]

    ## Calculate the expected and observed values per genotype
    ch=[]

    for g in exp[cross].keys():
        ## Get the expected number scaled to the number of non-missing genotypes
        e=exp[cross][g]*nonmis
        try:
            ## Calculate the deviation
            ch.append((obs[g]-e)**2/e)
        except KeyError:
            ## If an expected genotype was not observed
            ch.append((0-e)**2/e)


    ## Sum the genotype deviations to get chi square for the cross. 
    chi=sum(ch)

    ## Calculate the number of degrees of freedom, i.e. the number of expected genotypes
    ## possible in a given cross - 1. 
    n=len(exp[cross].keys())-1

    # Find the critical value for 95% confidence*
    m = 1 - (p/n_snps)
    crit = stats.chi2.ppf(q = m, df=n)

    if chi>crit:
        outcome=0   # reject
    else:
        outcome=1   # accept
        
    return outcome


def run_Mendel_family(p,gtp, parent1, parent2):
    ## Drop unnecessary columns:
    gtp = gtp.drop(columns=['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
    gtp = gtp.set_index('IID')
    
    ## Transpose the table - SNPs in rows, individuals in cols. 
    gtp_T = gtp.T
    start_SNPs = gtp_T.shape[0]
    
    ## Remove SNPs missing in either of the parents
    dfm = gtp_T.dropna(subset=[parent1,parent2])
    nm_snps=dfm.shape[0]
    
    ## Remove Mendelian inconsistencies
    dfm['Mendel']=dfm.apply(Mendel, args=(p,nm_snps, parent1, parent2), axis=1)
    
    ## Extract genotypes which passed the Mendelian inconsistencies check
    M = dfm.loc[dfm['Mendel']==1]
    ## Extract the list of rejected genotypes
    rejected=list(dfm.loc[dfm['Mendel']==0].index)
    
    ## Drop the Mendel column
    M=M.drop(columns=['Mendel'])
    
    ## Count number of SNPs retained
    mend_snps= M.shape[0]
    
    ## Count SNPs lost due to missingness
    missInParents = start_SNPs-nm_snps
      
#     print("Finished filtering and formating the genotypes.")
#     print("Starting number of SNPs: ", start_SNPs)
#     print("SNPs removed due to missing parental genotype: ", start_SNPs-nm_snps)
#     print("SNPs removed due to Mendelian inconsistencies: ", nm_snps-mend_snps)
#     print("Retained SNPs: {}".format(M.shape[0]))
    
    return M, rejected, missInParents

    
################################################################################################################################################
################################################################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument("--p", default=0.05, help="Specify p-value to be used for chi2 test, default 0.05")
parser.add_argument("--gtp", default="plink.raw", help="Specify the name of the file with genotypes, default plink.raw")
parser.add_argument("--ped", default="ped.csv", help="Specify the name of the pedigree file, default ped.csv")

args = parser.parse_args()

p=float(vars(args)['p'])
inFile=vars(args)['gtp']
ped=vars(args)['ped']

     
print("Input file: ", inFile) 
print("Pedigree file: ", ped)   
print("Probability threshold for chi-square test: ", p)

## Read in the pedigree and check how many families there are
ped = pd.read_csv(ped)

## Remove samples with unknown pedigree
ped = ped.dropna(subset=['Pedigree'])
fams = ped['Pedigree'].unique()

## Read in the genotypes
gtp = pd.read_table(inFile, delim_whitespace=True)

## If one family, run the usual
if len(fams)==1:
    ## Get parent IDs
    parent1 = ped['Mother.ID'].unique()[0]
    parent2 = ped['Father.ID'].unique()[0]

    ## Run the function on family
    rejected={}
    missP = {}
    M, rej, missInPar = run_Mendel_family(p, gtp, parent1, parent2)
    rejected["{}_{}".format(parent1, parent2)]=rej
    missP["{}_{}".format(parent1, parent2)]=missInPar
else:
    ## If more than one family
    ## Create an empty dataframe to store results - all SNPs included 
    snps = gtp.columns[6:]
    M = pd.DataFrame(index=snps)
    ## Create a dictionary where key will be family, value will be list of SNPs rejected
    rejects={}
    missP = {}
    ## For each family extract the parents and offspring
    for f in fams:
        print("\n")
        print("Running family {}".format(f))
        parent1 = ped.loc[ped['Pedigree']==f]['Mother.ID'].str.lower().unique()[0]
        parent2 = ped.loc[ped['Pedigree']==f]['Father.ID'].str.lower().unique()[0]
        family = list(ped.loc[ped['Pedigree']==f]['Sample.ID'])
        family.append(parent1)
        family.append(parent2)

        ## Extract the genotypes for the family
        gtpf = gtp.loc[gtp['IID'].isin(family)]

        ## Run the function within family
        try:
            Mtemp, rej, missInPar = run_Mendel_family(p, gtpf, parent1, parent2)
            ## Merge the results with the master frame
            test=M.merge(Mtemp, how='outer', left_index=True, right_index=True)
            M = test

            ## Save a list of rejected SNPs to dictionary
            rejects[f]=rej
            missP[f]=missInPar
            
        except KeyError:
            
            print("Failed to find parents for family ", f)
            print(parent1, parent2)

## Remove SNPs which failed in all families
M = M.dropna(axis=0, how="all")

## Get the number of retained SNPs
retained = M.shape[0]

## Format the genotypes to a raw format again
## Replace missing values with -9
M = M.fillna(-9)

## Change the format from float to int
M = M.astype('int32')

## Recode to raw format
t = M.T
t.reset_index(inplace=True)
t.insert(0,"FID", t['IID'])
t.insert(2,'PAT',0)
t.insert(3,'MAT',0)
t.insert(4,'SEX',0)
t.insert(5,'PHENOTYPE',-9)

## Save the results:
outgtp = inFile.replace(".raw","_MendelFiltered.ped")
outrej = inFile.replace(".raw","_MendelFiltered_rejectedSNPs.txt")
outmiss = inFile.replace(".raw","_MendelFiltered_SNPsMissingInParents.txt")

M.to_csv(outgtp, index=False, sep=' ')

r = open(outrej, 'w')
r.write("family\trejected_SNP\n")

for k,v in rejected.items():
    for snp in range(len(v)): 
        r.write("{}\t{}\n".format(k, v[snp]))
                
r.close()

print("Finished the Mendelian check for file: ", inFile)
print("{} families found".format(len(fams)))
print("Started with {} SNPs".format(len(snps)))
print("Removed {} SNPs as either missing in all parents, or failed in all families".format(len(snps)-M.shape[0])
print("Results saved to:")
print("