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
import pickle
import pandas as pd
import numpy as np
from itertools import chain
import time
#import cProfile
#from TEProb import transition_prob, emission_prob
from MetaMinimac import emission_prob, transition_prob, calcNumFlips, calcMetaDosages
from PosteriorProb import fwd_bwd
#from calcMetaDosages import calcMetaDosages
from calcDistMat import extract_int, calcLambda#, calcNumFlips
from IO import write_vcf, ds_gt_map
#from processData import sample_map
#Hidden states
#Hidden = (((1,1), (1,2)), ((2,1), (2,2)), ((1,1), (2,1)), ((1,1), (2,2)), ((1,2),(2,1)), ((1,2),(2,2)))
#Afr1, Afr2... Eur1, Eur2.. Afr1, Eur1..Afr1, Eur2..Afr2, Eur1..Afr2, Eur2, Afr2
#Hidden = (((1,1), (1,2)), ((2,1), (2,2)))
Hidden = (1, 2)
def sample_map(sampleID):
    return(dicto[sampleID]) #- 1) #index starts at 0


# %%
# #%run loadData.ipynb
#GL data path to file.. 
 
GL_path = "/net/fantasia/home/kiranhk/HMM/hapEASGL.csv"
#"ASWHaplotype0.csv"
#"ASWGL.csv" 30x
#"ASWGT_testcase.csv"
#"/net/fantasia/home/kiranhk/HMM/GL.csv"

#Allelic Dosages path to file...
  
AD_path = "/net/fantasia/home/kiranhk/HMM/230923_haploid_EASallelicdosages.npy"

#"230820_ASWallelicdosages_testcase.npy"
#"230813_ASWallelicdosages.npy" #ASW 30x 
#"/net/fantasia/home/kiranhk/HMM/230721_allelicdosages.npy" #aDNA


SNP_path = "230813ASWSNPs.npy"

#"/net/fantasia/home/kiranhk/HMM/230726SNPs.npy"

ad = np.load(AD_path)
gl = pd.read_pickle(GL_path)
#dist_mat = np.load(dist_path)
#remove semi-colons
SNPs = 'chr' + gl.index
#SNPs = np.load(SNP_path, allow_pickle = True)
dicto = pickle.load(open('EASdicto.p', 'rb'))
#dicto = {'I0054':1, 'I0103' :2, 'I0871':3, 'I0873':4}
assert ad.size/(2*2*len(dicto)) == len(gl) == len(SNPs) #check file size is consistent indicating markers are lined up
assert len(np.unique(SNPs))==len(SNPs) #check SNPs are unique

# %%
M=len(SNPs); L=30000
chunks = [np.arange(0+L*k,min(L*(k+1) + 1, M + 1)) for k in range(0, M//L + 1)]
print("Number of Chunks is ..", len(chunks))
#chunks

# %%
start = time.time() #start timing

samples = {}
#weights = {}

lda = calcNumFlips(calcLambda(SNPs)) #do this once and then subset

for sample in dicto.keys(): 
    mdosages = []
    #weightsc = []
    for c in chunks:
    #c = chunks[0]
        print("Meta Imputing sample ...", sample, "from", gl.iloc[min(c),:].name, "to", gl.iloc[max(c) - 1,:].name)
        print("Chunk size is...", max(c)-min(c))
    #subset data structures
        og_transformed = gl[min(c):max(c)][sample]
        #print(og_transformed)
        adc = ad[:, :, :, min(c):max(c)]
        ldac = lda[min(c):max(c),:]

    #calculate posteriors
    # #%timeit 
        pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)
    #cProfile.run('pst = fwd_bwd(Hidden, og_transformed.size, list(og_transformed), transition_prob, emission_prob, sample_map(sample), adc, ldac)')
        #print(pst)
    #calculate meta dosages
        mdosages.append(calcMetaDosages(pst, sample_map(sample), adc))
        #weightsc.append(pst)
    #add to samples
    samples[sample] = list(chain.from_iterable(mdosages))
    #weights[sample] = list(chain.from_iterable(weightsc))
   
   # np.save("metaASW" + sample + "230812", list(chain.from_iterable(mdosages)))
    #np.save("benchmark080923I0054", list(chain.from_iterable(mdosages)))
print("writing out vcf...")
write_vcf(samples, SNPs, "EAStargetphasingsimresults")
end = time.time ()
print("total time is", end - start)
#pickle.dump(samples, open('phasedsamples230919.p', 'wb')) #must use pickle to perserve dict


# %% [raw]
# #np.equal(
# a = np.load("benchmark080923I0054.npy")
# b = np.load("metaI0054230731_30kchunk.npy")
# c = np.load("benchmark080223I0054.npy")
# print(np.allclose(a,b, atol = 0.5))
# #len(a)
# #len(b)
# print(np.where(np.isclose(a,b, atol = 0.1)==False))
# np.allclose(a,c)

# %% [raw]
# weights = pickle.load(open('wts.p', 'rb'))
# #samples["NA19625"][0:20], SNPs[0:20]


# %% [raw]
# wts = pd.read_csv("~/1kg30xASW/weights.txt", sep = '\t')
# wts.index = 'chr' + wts.pop('#[1]ID')
# wts.pop("Unnamed: 123")
# truth = wts.iloc[:,0].apply(lambda x: float(str.split(x, ",")[0]))

# %% [raw]
# truth = pickle.load(open('wts.p', 'rb'))
# truth = truth.iloc[:,0:2]
# truth.columns = 'Afr', 'Eur'
# truth = truth/2

# %% [raw]
# est = pd.DataFrame.from_dict(weights["NA19625"])

# %% [raw]
# diff = np.isclose(truth, est.iloc[:,0], atol = 0.01)

# %% [raw]
# np.corrcoef(truth, est.iloc[:,0])**2

# %% [raw]
# truth.iloc[9143:9146],est.iloc[9143:9146,0]

# %% [raw]
# est[~diff], truth[~diff]

# %% [raw]
# df = pd.DataFrame(weights["NA19625"])
# df.index = SNPs
# df.columns = 'Afr', 'Eur'
# df.describe()

# %% [raw]
# import matplotlib.pyplot as plt
# plt.plot(df.index, df['Afr'], label='Afr')
# plt.plot(df.index, df['Eur'], label='Eur')
#
# plt.xlabel('Index')
# plt.ylabel('Values')
# plt.title('Weights across Markers for Sample NA19703')
# plt.legend()
#
# plt.show()
