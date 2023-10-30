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
import subprocess
from io import StringIO
import pickle
import json
import numpy as np
min_value = 1e-10  # Define the minimum value for the elements in GLs same as GLIMPSE2
epsilon = 1e-5


# %%
def ds_gt_map(ds): 
    ds = round(ds)
    if ds==0: 
        return ((0,0))
    if ds==1: 
        return ((0,1))
    if ds==2: 
        return((1,1))
    else: 
        raise ValueError("DS cannot be other than 0,1,2")


# %%
def write_vcf(samples, SNPs, outname, haploid=False): 
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
    vcfh.add_line("##FPLOIDY=1")

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
            if haploid: 
                r.samples[sample]['GT'] = round(samples[sample][s])
            else: 
                r.samples[sample]['GT'] = ds_gt_map(samples[sample][s])
            #r.samples[sample]['GP'] = (1,0,0)

        # Write this record to our VCF file
        vcf.write(r)

    # Close the VCF file
    vcf.close()


# %%
def get_sample_names(vcf_file):
    """Extract sample names from a VCF/BCF file and return a dictionary."""
    # The command to list sample names
    cmd = ["bcftools", "query", "-l", vcf_file]

    # Execute the command and get the result
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check for errors
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}")
        return None

    # Split the result by lines and form the dictionary
    samples = result.stdout.strip().split("\n")
    sample_dict = {sample: i+1 for i, sample in enumerate(samples)}

    return sample_dict

# Test the function
#vcf_path = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_haploid_chr20_ligated.bcf"  # or .vcf
#dicto = get_sample_names(vcf_path)


# %%
def query_bcftools(vcf_file, dicto, query_field):
    """Run bcftools and return the output as a pandas DataFrame."""
    cmd = [
        "bcftools", "query",
        "-f", f"%CHROM:%POS:%REF:%ALT\t[%{query_field}\t]\n", 
        vcf_file
    ]

    # Execute the command and get the result
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    # Check for errors
    if result.returncode != 0:
        print(f"Error running bcftools: {result.stderr}")
        return None

    # Convert the result to a DataFrame
    df = pd.read_csv(StringIO(result.stdout), sep="\t", header=None)
    df.columns = ["ID"] + list(dicto) + ["drop"]

    return df
# Test the function
#vcf_path = "/net/fantasia/home/kiranhk/software/glimpse/tutorial/GLIMPSE_ligate/AFR_haploid_chr20_ligated.bcf"
#df = query_bcftools(vcf_path, dicto, "DS")
#print(df.head())


# %%
def unphred(GL_tuple):
    # Create your list comprehension, incorporating the minimum value check
    GLs = [max(pow(10., -i/10.), min_value) for i in GL_tuple]
    
    return tuple(GLs)  # Convert the list to a tuple before returning


# %%
def read_vcfs(GL_path, DS_path1, DS_path2): 
    #get PL --> unphred ---> pandas df #get DS ---> create allelic dosage structure numpy
    
 #file checks sample names are consistent
    assert get_sample_names(DS_path1) == get_sample_names(GL_path) == get_sample_names(DS_path2)
 
 #use bcftools to extract data
    dicto = get_sample_names(DS_path1)
    gl = query_bcftools(GL_path, dicto, "PL")
    
    #get ploidy
    ploidy = len(json.loads("[" + gl.iloc[1,1] + "]")) - 1 #assuming gl all have the same ploidy
    print("ploidy is ...", ploidy)
    

    
    if ploidy ==1:
        eur_dosage = query_bcftools(DS_path2, dicto, "DS")
        afr_dosage = query_bcftools(DS_path1, dicto, "DS")
    else: 
        eur_dosage = query_bcftools(DS_path2, dicto, "AP")
        afr_dosage = query_bcftools(DS_path1, dicto, "AP")
        
   

    #embedded function that knows what dicto is, otherwise its undefined
    def sample_map(sampleID):
        return(dicto[sampleID] - 1)     
        


    afr_dosage.pop("drop")
    eur_dosage.pop("drop")
    gl.pop("drop")

#merge Dosage files
    #all_dosage = pd.merge(eur_dosage, afr_dosage, on="ID")
#Need to consider all variants now: 
    all_dosage = pd.merge(eur_dosage, afr_dosage, on="ID", how='outer')
    
#OR 
 #inner join is imputed and outer - inner is saved to a vcf file which is then (merged?) to the one generated by the output function
#write a function for this and test it    
    
    
#left merge GL 
    all_GL = pd.merge(all_dosage["ID"], gl, how = 'left', on = "ID")
    all_GL.index = all_GL.pop("ID") #should really be merge GL
    all_dosage.index = all_dosage.pop("ID")
    #print(all_dosage)

#all dosage str --> tuple    
    if ploidy==2: #for 0 dosage case
        all_dosage = all_dosage.applymap(lambda x: (0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
    
#GL tuple, unphred
    if ploidy == 1: 
        all_GL = all_GL.applymap(lambda x: (0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
    elif ploidy == 2: 
        all_GL = all_GL.applymap(lambda x: (0,0,0) if pd.isna(x) else tuple(json.loads("[" + x + "]")))
    pos = (all_GL.applymap(len)).apply(any, axis = 1)
    assert all(pos) #check for len = 0 "truthy" means 0 is false and everything else is true
    assert np.all(all_GL.index == all_dosage.index) #check SNPs
    SNPs = all_dosage.index #output 1

    obsGLMerge=all_GL.applymap(unphred)

    print("obsGL:", len(gl), "left merged GL:", len(all_GL), "merged dosage:",len(all_dosage), "african:", len(afr_dosage), "eur:", len(eur_dosage))
    assert len(all_dosage) == len(all_GL) and len(afr_dosage) > len(eur_dosage) and len(afr_dosage) < len(all_dosage) and len(eur_dosage) < len(all_dosage) 

    #create allelic dosages
    M=all_dosage.shape[0] 
    N = len(dicto)

    edosage = all_dosage.iloc[:,0:N] #2
    adosage = all_dosage.iloc[:,N:N*2] #1

    assert sum(edosage.index==adosage.index)==M
    assert sum(edosage.index==all_GL.index)==M #check SNPs are the same in dosages and GLs

    edosage.columns = list(dicto.keys())
    adosage.columns = list(dicto.keys())


    sample_list = dicto.keys()

    allelic_dosages = np.zeros((2, 2, N, M)) #could change this to be smaller for ploidy 1 
    

#loop over samples
    for sample in sample_list:
        if ploidy ==1: 
            allelic_dosages[1][0][sample_map(sample)] = edosage[sample].apply(lambda x: x) #2
            allelic_dosages[0][0][sample_map(sample)] = adosage[sample].apply(lambda x: x) #1
        elif ploidy == 2: 
            allelic_dosages[1][0][sample_map(sample)] = edosage[sample].apply(lambda x: x[0]) #(2,1)
            allelic_dosages[0][0][sample_map(sample)] = adosage[sample].apply(lambda x: x[0]) #(1,1)
            allelic_dosages[1][1][sample_map(sample)] = edosage[sample].apply(lambda x: x[1])#(2,2)
            allelic_dosages[0][1][sample_map(sample)] = adosage[sample].apply(lambda x: x[1])#(1,2)

   #clip values to prevent underflow
    ad_clipped = np.clip(allelic_dosages, epsilon, 1 - epsilon)
            
    return(SNPs, dicto, obsGLMerge, ad_clipped)


# %% [raw]
# GL = "/net/fantasia/home/kiranhk/1kg30xASW/230810genogvcfs.vcf.gz"
#
# DS_afr = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_chr20_ligated.bcf"
#
# DS_eur = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUR_chr20_ligated.bcf"
#
# test = read_vcfs(GL, DS_afr, DS_eur)

# %% [raw]
# dicto_test = pickle.load(open("dicto.p", "rb"))
# print(test[1] == dicto_test)
#
# SNPs_test = np.load("230913ASWSNPs.npy", allow_pickle = True)
#
# print(np.all(test[0] == SNPs_test))
#
# #test allelic dosages and obs gl 
# #ad_test = np.load("230918_haploid_ASWallelicdosages.npy")
# #np.allclose(test[3], ad_test)

# %% [raw]
# gl_test = pickle.load(open("ASWGL.csv", "rb"))
# #gl_test.index = test[2].index
#

# %% [raw]
# test[2].describe()

# %% [raw]
# import pickle
# import numpy as np
# samples = pickle.load(open('ASWmetaimputed_allsamples.p', 'rb'))
# SNP_path = "230813ASWSNPs.npy"
# SNPs = np.load(SNP_path, allow_pickle = True)

# %% [raw]
# write_vcf(samples, SNPs, "example")
