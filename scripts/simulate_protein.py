#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 12:12:01 2022

@author: ptruong
"""

import numpy as np

number_of_replicates = 3
number_of_conditions = 2 
number_of_shared_peptides = 4

number_of_uniq_peptides_per_proteoform = 4 # make a list

number_of_uniq_peptides_per_proteoform = [3, 3, 3]

prot1_per_condition = [1.0, 1.0] # A B - concentration per condition
prot2_per_condition = [0.3, 3.0]

proteoform_per_condition = [[1.0, 1.0],
                            [0.3, 3.0],
                            [1., 2.]]

ionization_mu = 3.
ionization_sigma = 5. 
noise_proportionality = 0.2
noise_min = 0.1


assert len(proteoform_per_condition) == len(number_of_uniq_peptides_per_proteoform), f"len(proteoform_per_condition) != len(proteoform_per_condition) i.e. ({len(proteoform_per_condition)} != {len(number_of_uniq_peptides_per_proteoform)})"



number_of_peptides = number_of_uniq_peptides_per_proteoform*2 + number_of_shared_peptides
rng = np.random.default_rng()    

def add_noise(raw):
    stdv = max(noise_min, raw * noise_proportionality)
    noise = rng.normal(0., stdv)
    return abs(raw + noise) # guard ourseves from negative values
    
add_noise_v = np.vectorize(add_noise)

[[proteoform_per_condition[i]*number_of_uniq_peptides_per_proteoform + \
  proteoform_per_condition[i]*] for i in range(len(proteoform_per_condition))]

n_proteoforms = len(proteoform_per_condition)
n_conditions = len(proteoform_per_condition[0])

peptide_conc_per_condition = []
for i in range(n_conditions):
    condition_peptide_conc = []
    for j in range(n_proteoforms):
        condition_peptide_conc += [proteoform_per_condition[j][i]]*number_of_uniq_peptides_per_proteoform[j]
    peptide_conc_per_condition.append(condition_peptide_conc)
            
    
# Add shared peptides to peptide_conc_per_condition
[prot1_per_condition[0]+prot2_per_condition[0]]*number_of_shared_peptides


peptide_conc_A = np.asarray([prot1_per_condition[0]]*number_of_uniq_peptides_per_proteoform + \
                 [prot1_per_condition[0]+prot2_per_condition[0]]*number_of_shared_peptides + \
                 [prot2_per_condition[0]]*number_of_uniq_peptides_per_proteoform)

peptide_conc_B = np.asarray([prot1_per_condition[1]]*number_of_uniq_peptides_per_proteoform + \
                 [prot1_per_condition[1]+prot2_per_condition[1]]*number_of_shared_peptides + \
                 [prot2_per_condition[1]]*number_of_uniq_peptides_per_proteoform)

peptide_conc = np.vstack((np.tile(peptide_conc_A, (number_of_replicates, 1)),\
                          np.tile(peptide_conc_B, (number_of_replicates, 1))))

ionization = rng.lognormal(ionization_mu, ionization_sigma, number_of_peptides)
peptide_abundance_raw = np.multiply(ionization,peptide_conc)
peptide_abundance = add_noise_v(peptide_abundance_raw)    








    