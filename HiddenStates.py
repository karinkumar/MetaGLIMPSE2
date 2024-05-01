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
def generate_hidden(K, mixed, haploid): 
    '''
    K is the number of reference panels
    mixed states is a boolean signifying whether mixed states are allowed
    haploid is a boolean signifying ploidy
'''
    if mixed and haploid: 
        raise ValueError("Cannot have states that are both mixed and haploid")
    tuples = []
    if haploid: #cannot be both haploid and mixed 
        return tuple(range(1, K+1)) 
    if not haploid: 
        for k in range(1, K+1): 
            tuples.append(((k, 1), (k, 2))) #handle case where reference panel is the same
        if not mixed:    
            return(tuple(tuples))
        else: 
            for i in range(1, K):
                for j in range(i + 1, K + 1):
                    tuples.append(((i, 1), (j, 1)))
                    tuples.append(((i, 1), (j, 2)))
                    tuples.append(((i, 2), (j, 1)))
                    tuples.append(((i, 2), (j, 2)))

        return tuple(tuples)        

# %% [raw]
# #Hidden states the way it was
# #if not mixed_states and not haploid: 
# Hidden_plain = (((1,1), (1,2)), ((2,1), (2,2)))
# #elif haploid and not mixed_states: 
# Hidden_hap = (1, 2)
# #else: 
# Hidden_mix = (((1,1), (1,2)), ((2,1), (2,2)), ((1,1), (2,1)), ((1,1), (2,2)), ((1,2),(2,1)), ((1,2),(2,2)))
#     
# #Afr1, Afr2... Eur1, Eur2.. Afr1, Eur1..Afr1, Eur2..Afr2, Eur1..Afr2, Eur2, Afr2


# %% [raw]
# #test cases
# assert generate_hidden(2, True, False) == Hidden_mix
# assert generate_hidden(2, False, False) == Hidden_plain
# assert generate_hidden(2, False, True) == Hidden_hap
