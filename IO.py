# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import pysam
from pysam import VariantFile
import pandas as pd


# %%
def ds_gt_map(ds): 
    ds = round(ds)
    if ds==0: 
        return((0,0))
    if ds==1: 
        return((0,1))
    if ds==2: 
        return((1,1))
    else: 
        raise ValueError("DS cannot be other than 0,1,2")


# %%
def write_vcf(samples, SNPs, outname): 
#write vcf
    vcfh = pysam.VariantHeader()
# Add a sample named "ahstram" to our VCF header
    for s in samples.keys():
        vcfh.add_sample(s)
# Add a contig (chromosome 20) to our VCF header
    vcfh.add_meta('contig', items=[('ID','chr20')])
    
# Add GT andDS to FORMAT in our VCF header
    vcfh.add_meta('FORMAT', items=[('ID',"GT"), ('Number',1), ('Type','String'),
    ('Description','Best Guess Genotype')])
    vcfh.add_meta('FORMAT', items=[('ID',"DS"), ('Number','A'), ('Type','Float'),
    ('Description','Meta Imputed Genotype Dosage')])
    vcfh.add_line("##FPLOIDY=2")

    vcfh.add_meta('FORMAT', items=[('ID',"GP"), ('Number',3), ('Type','Float'),
    ('Description','Meta Imputed Genotype Dosage')])
    
# Open a file, "example.vcf" for writing, with our "vcfh" header
    vcf = pysam.VariantFile(outname + ".vcf.gz", "w", header=vcfh)

# Create a record at chr1:1000 A/T which failed filter due to "RF"
# The 'start' value is 0-based, 'stop' is 1-based
    for s,val in enumerate(SNPs): 
        ID = str.split(val, ":")
        r=vcf.new_record(contig='chr20', start=int(ID[1]) - 1 , stop=int(ID[1]),
            alleles=[ID[2],ID[3]])
        for sample in samples.keys():
        # Set dosage
            r.samples[sample]['DS'] = samples[sample][s]
            r.samples[sample]['GT'] = ds_gt_map(samples[sample][s])
            r.samples[sample]['GP'] = (1,0,0)

        # Write this record to our VCF file
        vcf.write(r)

    # Close the VCF file
    vcf.close()


# %%
def read_vcfs(GL_path, DS_path1, DS_path2): 
    #get PL --> unphred ---> pandas df
    #get DS ---> create allelic dosage structure numpy
GL = ...
DS_afr = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_haploid_chr20_ligated.bcf"
DS_eur = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUR_haploid_chr20_ligated.bcf"
    
    
#use bcftools to extract data

getPL(GL)
getDS(DS)

#get sample names via bcftools
assert getNames(DS) == getNames(GL) == getNames(DS2)
dicto = getNames(GL)

#merge Dosage files

#left merge GL 

#GL tuple, unphred

#DS tuple


#create allelic dosages

#file checks

#export files(how?, return them as a list?)

# %%
def getPL(): 
...
def getDS():
...
def getNames(): 

# %% [raw]
# import pickle
# import numpy as np
# samples = pickle.load(open('ASWmetaimputed_allsamples.p', 'rb'))
# SNP_path = "230813ASWSNPs.npy"
# SNPs = np.load(SNP_path, allow_pickle = True)

# %% [raw]
# write_vcf(samples, SNPs, "example")
