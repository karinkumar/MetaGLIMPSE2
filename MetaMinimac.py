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
import numpy as np
H=2


# %%
def calcNumFlips(lda): 
    arr = np.zeros((lda.size, 3))
    arr[:,2] = lda/2
    arr[:,0] = 1 - lda/2 #numflips is 0
    return(arr)


# %%
def emission_prob(hidden, obs, m, sampleID, npa):
    '''Hidden is Afr 1, Eur 2. Always using 0th haplotype. 
    '''
    #print(hidden)
    d = npa[hidden - 1][0][sampleID][m] + 1e-5
   
    GL0, GL1 = obs
    #print("hidden", hidden, "dosage", d, "GL0", GL0, "GL1", GL1, "emission", ((1 - d)**2)*GL0 + (d**2)*GL1)
    return(((1 - d)**2)*GL0 + (d**2)*GL1) #square to avoid non-informative prior
    
    #if obs==0: 
        #print("allele count", obs, "hidden", hidden, "marker", m, "emissions",  1 - d)
      #  return(1 - d + 1e-5)
        
    #elif obs==1: 
        #print("allele count", obs, "hidden", hidden, "marker", m, "emissions",  d)
     #   return d + 1e-5
        
    #else: 
        #raise ValueError("A must be 0 or 1")


# %%
def transition_prob(from_state, to_state, lda, m): 
    '''
    Input the to/fro hidden states (each a tuple) e.g. ((1,1),(1,2)) is Ref1Allel1,Ref1Allele2 the distance matrix, and marker (i.e time T along the HMM)
    Ouput transition probability
   '''
    num_flips = 2 - 2*int(from_state == to_state)
    #print("from state", from_state, "to state", to_state, "numflips", num_flips, "lda", lda[m, num_flips], "marker", m)
    return(lda[m, num_flips])


# %%
def calcMetaDosages(posteriors, sample, allelic_dosages): 
    meta_dosages=list()
    for m in range(len(posteriors)): 
        panels = posteriors[m]
        meta_dosage=0
        for key, value in panels.items():
            a = key #key[0]
            #c = key[1]
            meta_dosage += allelic_dosages[a-1][0][sample][m]*value #+ allelic_dosages[c-1][0][sample][m]*value
        meta_dosages.append(meta_dosage)
        meta_dosage=0 #reset weighted sum
    #print(max(meta_dosages), min(meta_dosages))
    if min(np.round(np.array(meta_dosages),3))>= 0 and max(np.round(np.array(meta_dosages),3)) <= 1: 
        return(meta_dosages)
    else: 
        raise ValueError("Meta Dosages Not Between 0 and 1")

# %% [raw]
# import pickle
# import pandas as pd
# import numpy as np
# #read in posterios 
# weights = pickle.load(open('allphasedposteriors.p', 'rb'))
# #pick 10 posteriors at random, hand calculate dosages
# weights
# est = pd.DataFrame.from_dict(weights)
# #read in allelic dosages
# ad = np.load("/net/fantasia/home/kiranhk/HMM/230918_haploid_ASWallelicdosages.npy")
# #snps to get ID for checking

# %% [raw]
# est.describe()

# %% [raw]
# for m in range(280929): 
#     print(ad[0,0,60,m]*0.5 + ad[1,0,60,m]*0.5)
