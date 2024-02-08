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
###command line interface###

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--afr", required=True)
parser.add_argument("--eur", required = True)
parser.add_argument("--gl", required = True)
parser.add_argument("--out", required = True)
parser.add_argument("--haploid", required = False, dest = 'haploid', action = 'store_true')
parser.add_argument("--nomixedstate", required = False, dest = 'nomixedstate', action = 'store_true')
parser.add_argument("--pickle", required = False, dest = 'pickle', action = 'store_true')
args = parser.parse_args()
parser.set_defaults(haploid=False)
parser.set_defaults(nomixedstate=False)
parser.set_defaults(pickle=False)
K=2 #count number of args
haploid = args.haploid #False
mixed_states = not args.nomixedstate #False
if haploid and mixed_states: #sanity check override user 
    raise ValueError("Cannot have mixed states for haploid data")
print("mixed states are...", mixed_states)

GL = args.gl #"/net/fantasia/home/kiranhk/1kg30xEAS/genogvcfs1x.vcf.gz"

DS_afr = args.afr #"/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/AFR_EASdiploid_chr20_ligated.bcf"

DS_eur = args.eur #"/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUR_EASdiploid_chr20_ligated.bcf"

# %% [raw]
#

# %% [raw]
# ###run through notebook interface###
# K=2 #number of hidden states
# haploid = False
# mixed_states = False
# if haploid and mixed_states: #sanity check override user 
#     raise ValueError("Cannot have mixed states for haploid data")
# print("mixed states are...", mixed_states, "... with", K, "reference panels")
#
# GL ="/net/fantasia/home/kiranhk/1kg30xEUR/gl/bcftoolsgenogvcfs6x.vcf.gz"
#
# DS_afr = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUREURAdiploid_6xchr20.vcf.gz"
#
# DS_eur = "/net/fantasia/home/kiranhk/software/GLIMPSE2_for_kiran_kumar/GLIMPSE_ligate/EUREURBdiploid_6xchr20.vcf.gz"

# %%
import pickle
import pandas as pd
import numpy as np
from itertools import chain
import time
#import cProfile
from PosteriorProb import fwd_bwd
if haploid:
    from MetaMinimac import emission_prob, transition_prob, calcNumFlips, calcMetaDosages
else: 
    from TEProb import transition_prob, emission_prob
    from calcMetaDosages import calcMetaDosages
if haploid: 
    from calcDistMat import extract_int, calcLambda
else: 
    from calcDistMat import extract_int, calcLambda, calcNumFlips
from IO import write_vcf, ds_gt_map, read_vcfs
from HiddenStates import generate_hidden
Hidden = generate_hidden(K, mixed_states, haploid)    
def sample_map(sampleID):
    return(dicto[sampleID] - 1) #index starts at 0


# %%
print("Reading vcfs ...")
 #start timing
start = time.time()
SNPs, dicto, gl, ad, all_dosage = read_vcfs(GL, DS_afr, DS_eur)

# %%
print("Checking vcfs...")
assert ad.size/(2*2*len(dicto)) == len(gl) == len(SNPs) #check file size is consistent indicating markers are lined up
assert len(np.unique(SNPs))==len(SNPs) #check SNPs are unique
assert len(dicto) == gl.shape[1] #check sample names are equivalent

M=len(SNPs); L=30000
chunks = [np.arange(0+L*k,min(L*(k+1) + 1, M + 1)) for k in range(0, M//L + 1)]
print("Number of Chunks is ..", len(chunks))
#chunks

# %%
samples = {}
#weights = {}

lda = calcNumFlips(calcLambda(SNPs)) #do this once and then subset

for sample in dicto.keys(): 
    mdosages = []
    #weightsc = []
    for c in chunks:
        print("Meta Imputing sample ...", sample, "from", gl.iloc[min(c),:].name, "to", gl.iloc[max(c) - 1,:].name)
        print("Chunk size is...", max(c)-min(c))
    #subset data structures
        og_transformed = gl[min(c):max(c)][sample]
        adc = ad[:, :, :, min(c):max(c)]
        ldac = lda[min(c):max(c),:]
 

    #calculate posteriors
        # #%timeit 
        pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)
    #cProfile.run('pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)')
    #calculate meta dosages
        #np.save("231016posteriors_test", pst)
        mdosages.append(calcMetaDosages(pst, sample_map(sample), adc))
        #weightsc.append(pst)
    #add to samples
    samples[sample] = list(chain.from_iterable(mdosages))
    #weights[sample] = list(chain.from_iterable(weightsc))
   


# %% [raw]
# #plot weights
# df = pd.DataFrame(weights["HG00121"])
# #df.index = SNPs
# df.columns = 'EURA', 'EURB'
# df.describe()
#
# import matplotlib.pyplot as plt
# plt.plot(df.index, df['EURA'], label='EURA')
# plt.plot(df.index, df['EURB'], label='EURB')
#
# plt.xlabel('Index')
# plt.ylabel('Values')
# plt.title('Weights across Markers for Eur Target Sample HG00121')
# plt.legend()
#
# plt.show()

# %%
    
print("writing out vcf...")
write_vcf(samples, SNPs, args.out)
end = time.time ()
print("total time is", end - start)
if args.pickle:
    pickle.dump(samples, open(args.out + '.p', 'wb')) #must use pickle to perserve dict


# %% [raw]
# %run loadData.ipynb
# #GL data path to file.. 
#  
# GL_path = "/net/fantasia/home/kiranhk/HMM/hapASWGL.csv"
# #"ASWHaplotype0.csv"
# #"ASWGL.csv" 30x
# #"ASWGT_testcase.csv"
# #"/net/fantasia/home/kiranhk/HMM/GL.csv"
#
# #Allelic Dosages path to file...
#   
# AD_path = "/net/fantasia/home/kiranhk/HMM/230918_haploid_ASWallelicdosages.npy"
#
# #230923_haploid_EASallelicdosages.npy" #EAS 30x
# #"230820_ASWallelicdosages_testcase.npy"
# #"230813_ASWallelicdosages.npy" #ASW 30x 
# #"/net/fantasia/home/kiranhk/HMM/230721_allelicdosages.npy" #aDNA
#
#
# SNP_path = "230913ASWSNPs.npy" #also the same for EAS
#
# #"/net/fantasia/home/kiranhk/HMM/230726SNPs.npy"
#
# ad = np.load(AD_path)
# gl = pd.read_pickle(GL_path)
# #dist_mat = np.load(dist_path)
# #remove semi-colons
# #SNPs = 'chr' + gl.index
# SNPs = np.load(SNP_path, allow_pickle = True)
# dicto = pickle.load(open('dicto.p', 'rb'))
