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

def calcMetaDosages(posteriors, sample, allelic_dosages):
    '''INPUT: posteriors list of dictionaries, sample number 0 - N, allelic dosages numpy array
       OUTPUT: list of meta genotype dosages one per marker 0 to M
       
       NOTES:multiply across haplotypes that have the SAME posterior weight 1,1 represents reference panel 1 allele 1  
       AND reference panel 2 allele 2 (d_a*0.01 +  d_b*0.01)
    '''
    meta_dosages=list()
    for m in range(len(posteriors)): 
        panels = posteriors[m]
        #print(min(panels), max(panels))
        meta_dosage=0
        for key, value in panels.items():
            a,b = key[0]
            c,d = key[1]
            meta_dosage += allelic_dosages[a-1][b-1][sample][m]*value + allelic_dosages[c-1][d-1][sample][m]*value
        meta_dosages.append(meta_dosage)
        meta_dosage=0 #reset weighted sum
    #print(max(meta_dosages), min(meta_dosages))
    if min(np.round(np.array(meta_dosages),3))>= 0 and max(np.round(np.array(meta_dosages),3)) <= 2: 
        return(meta_dosages)
    else: 
        print(min(np.round(np.array(meta_dosages),3)), max(np.round(np.array(meta_dosages),3)) )
        raise ValueError("Meta Dosages Not Between 0 and 2")

# %% [raw]
# posteriors = np.load("231016posteriors_test.npy", allow_pickle = True)
# allelic_dosages = np.load("230813_ASWallelicdosages.npy", allow_pickle = True)

# %% [raw]
# sample = 8
#
# meta_dosages=list()
# for m in [26499]: #range(len(pst)): 
#     panels = posteriors[m]
#     print(panels)
#     meta_dosage=0
#     for key, value in panels.items():
#         a,b = key[0]
#         c,d = key[1]
#         meta_dosage += allelic_dosages[a-1][b-1][sample][m]*value + allelic_dosages[c-1][d-1][sample][m]*value
#     meta_dosages.append(meta_dosage)
#     meta_dosage=0 

# %% [raw]
# np.where(np.array(meta_dosages) < 0)
