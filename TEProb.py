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
import pandas as pd
#epsilon = 1e-5

# %% [raw]
# #Hidden states
# Hidden = (((1,1), (1,2)), ((2,1), (2,2)), ((1,1), (2,1)), ((1,1), (2,2)), ((1,2),(2,1)), ((1,2),(2,2)))
# #Afr1, Afr2... Eur1, Eur2.. Afr1, Eur1..Afr1, Eur2..Afr2, Eur1..Afr2, Eur2, Afr2
#
# num_flips_tbl = np.zeros((len(Hidden), len(Hidden))) - 8
# for num, state in enumerate(Hidden):
#    # print(num == Hidden.index(state))
#     a,b = state
#     for num_, state_ in enumerate(Hidden): 
#         c,d = state_ #matrix of from to pairs --> num flips
#         num_flips_tbl[num, num_] = 2 - (int(a==c or a==d) + int(b==d or b==c))
# num_flips_tbl
# nfpd = pd.DataFrame(num_flips_tbl, index=Hidden, columns=Hidden)

# %% [raw]
# #test num_flips regardless of order this should work
#
# from_state = ((1,1),(1,2))
# a,b = from_state
# a
# b
# to_state = ((1,2), (1,1))
# c,d = to_state
#
# 2 - (int(a==c or a==d) + int(b==d or b==c))
#

# %%
def transition_prob(from_state, to_state, lda, m): 
    '''
    Input the to/fro hidden states (each a tuple) e.g. ((1,1),(1,2)) is Ref1Allel1,Ref1Allele2 the distance matrix, and marker (i.e time T along the HMM)
    Ouput transition probability
   '''
    a,b = from_state
    c,d = to_state
    num_flips = 2 - (int(a==c or a==d) + int(b==d or b==c))
    #print(num_flips)
    #print(lda[m, num_flips])
    return(lda[m, num_flips])
    


# %% [raw]
# if num_flips ==0:
#         return (1 - lda)**2 + (2*(lda - lda**2))/H + lda**2/(H**2)
#     elif num_flips==1:
#         return ((1 - lda)*lda)/H + (lda**2)/H**2
#     elif num_flips==2:
#         return (lda)/(H)
#     else:
#         raise ValueError("Cannot have more than 2 flips")

# %% [raw]
# lda = np.load("ldatestcase.npy")
# #transition_prob(((1,1), (2,2)), ((1,2), (2,1)), lda, 1) + transition_prob(((1,1), (2,2)), ((1,1), (2,2)), lda, 1)+ 4*transition_prob(((1,1), (2,2)), ((1,1), (2,1)), lda, 1)   
# lda[0], lda[1], lda[2]
# #np.where(np.isclose(np.sum(lda, axis = 1),1, atol = 0.01)==False)
# #np.sum(lda[9144])

# %% [raw]
# unphred(float('inf'))
# #test (off by a constant so it doesn't matter?)
# unphred(100)/unphred(90), unphred(25)/unphred(15), unphred(10)/unphred(0)

# %%
def emission_prob(hidden, obs, m, sampleID, npa):
    '''INPUT: dosage numpy array (npa), observed genotype: 0, 1, 2, 3(missing), hidden state = ((1,1),(1,1)), sampleID, dosages data frame
        RETURNS: emission probability 
    '''
    #make sure obs is a tuple: #computationallly expensive? only use for debugging
    #if not isinstance(obs, tuple): 
     #   raise TypeError("Obs not Tuple at", m)
    #elif len(obs)!= 3:
     #   raise ValueError("Obs must be Tuple of Length 3 at", m)

    GL0, GL1, GL2 = obs #unpack
   # print(GL1)
    #No only missing genotype likleihood. We are only combining inputation sets.
    #Nothing is "missing anymore" each marker has a dosage. Also GL cannot be "missing" just has a flat prior
    a,b = hidden
    #print(a[0]-1, a[1]-1, sampleID, m)
    d_a = npa[a[0] - 1][a[1] - 1][sampleID][m] #dosage for hidden ref haplotype a 
    d_b = npa[b[0] - 1][b[1] - 1][sampleID][m]  #dosage for hidden ref haplotype b
  
   # d_a = np.clip(d_a, epsilon, 1 - epsilon) #clip to ensure no underflow issues
   # d_b = np.clip(d_b, epsilon, 1 - epsilon)
    
    return GL2 * d_a*d_b + GL1 *(d_a * (1 - d_b) + d_b * (1 - d_a)) + GL0* (1 - d_a) * (1 - d_b)

# %% [raw]
# ad = np.load("230820_ASWallelicdosages_testcase.npy")
# import pickle
# gl = pickle.load(open("ASWGT_testcase.csv", 'rb'))
# gl = gl["NA19625"]

# %% [raw]
# #dist_mat = np.load("dist_matchunk1.npy")
# #ad = np.load("/net/fantasia/home/kiranhk/HMM/230721_allelicdosages.npy")
# #
#
# #emission_prob(((1,1), (1,2)), gl[0], 0, 0, ad) #== 0.978001
# #emission_prob(((2,1), (2,2)), gl[0], 0, 0, ad)
# #emission_prob(((1,1), (1,2)), gl[1], 1, 0, ad)
# #emission_prob(((2,1), (2,2)), gl[1], 1, 0, ad)
# #emission_prob(((1,1), (1,2)), gl[2], 2, 0, ad)
# Hidden = (((1,1), (1,2)), ((2,1), (2,2)), ((1,1), (2,1)), ((1,1), (2,2)), ((1,2),(2,1)), ((1,2),(2,2)))
# for h in Hidden:
#     print(h, emission_prob(h, gl[1997], 1997, 0, ad)) #all tests fine
# #import time
# %timeit emission_prob(((1,1), (2,2)), (1,1,1), 1, 1, ad)
#
# %load_ext line_profiler 2.41e-06 s
# %lprun -f transition_prob transition_prob(((1,1), (2,2)), ((1,1), (2,2)), dist_mat, 1)
# %lprun -f emission_prob emission_prob(((1,1), (2,2)), (1,1,1), 1, 1, ad)
