# for comparison to Run_MSPrime_SingleRep_Dros.py, this rescales by reducing population size and split times by 50x while increasing mutation and recombination/conversation rates 50x
# three population split A,(B,C)

import concurrent.futures
import msprime
import sys
import os
import numpy as np
import math
import statistics as stats
import random

# Fst estimate follows Reynolds et al 1983 (Genetics 105:767) but for haploid rather than diploid samples
# i.e. n1,n2 replaced by n1/2, n2/2 etc in formulas from 769 in Reynolds et al
# some terms rearranged for simplicity, i.e. Numa as (p1-p2)^2 rather than (1/2)((p1-p2)^2 + (p2-p1)^2 etc
def rey_fst(freq1, freq2, SampleSize1, SampleSize2):
    n1 = SampleSize1
    n2 = SampleSize2
    #sec_allele_1 = 1 - freq1
    #sec_allele_2 = 1 - freq2
    SharedNum = n1*(freq1 - freq1**2) + n2*(freq2 - freq2**2)
    # check if a factor of 1/2 is needed for NumA, as per Reynolds et al 1983
    NumA = (freq1 - freq2)**2
    FracNum = ((n1 + n2)/2)*SharedNum
    FracDen = n1*n2*((n1+n2)/2 - 1)
    frac = FracNum/FracDen
    WholeNum = NumA - frac
    DenFracNum = (n1*n2 - (n1+n2)/2)*SharedNum
    DenFrac = DenFracNum/FracDen
    WholeDen = NumA + DenFrac
    if WholeDen !=0:
        FST = WholeNum/WholeDen
    else:
        FST = 0
    #if (FST > 1):
        #FST = 1
    #elif (FST < 0):
        #FST = 0
    #else:
        #FST = FST
    return([FST, WholeNum, WholeDen])
    
#

# calculate numerator and denominator terms of Reynolds Fst for every locus          
def multilocus_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
    fst_list = [rey_fst(Freqlist1[i], Freqlist2[i], SampleSize1, SampleSize2) for i in range(len(Freqlist1))]
    return(fst_list)


#
    
# Window Fst
#weighted average of multilocus fst: sum of numerator and denominator across sites    
def weighted_avg_rey_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
    all_dat = multilocus_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2)
    al_list = [all_dat[i][1] for i in range(len(Freqlist1))]
    albl_list = [all_dat[i][2] for i in range(len(Freqlist1))]
    s_albl = sum(albl_list)
    s_al = sum(al_list)
    if(s_albl > 0):
        temp_ratio = s_al/s_albl
    else:
        temp_ratio = 0
    if temp_ratio < 0:
        temp_ratio = 0
    elif temp_ratio > 1:
        temp_ratio = 1
    else:
        temp_ratio = temp_ratio
    return(temp_ratio)
    
# PBS = Population Branch Statistic (Yi et al 2010)
# Fst1: AB, Fst2: AC, Fst3: BC
# T: rescale as -log(1-Fst) for approximate additivity
# PBS branch length estimate, using same methodology as NJ   
def PBS_Fun(Fst_1, Fst_2, Fst_3):
    T1 = -math.log(1-Fst_1)
    T2 = -math.log(1-Fst_2)
    T3 = -math.log(1-Fst_3)
    PBS = (T1 + T2 - T3)/2
    return([PBS,T1,T2,T3])


# Population branch estimate PBE - rescale PBS wrt median PBE and T. Yassin et al (2016)
# PBS1: A focal lineage, 2:B, 3:C etc
#def PBE_Fun(Fst_1, Fst_2, Fst_3):
#    #PBS1: A focal, PBS2: B focal, PBS3: C focal
#    PBS1 = PBS_Fun(Fst_1, Fst_2, Fst_3)
#    PBS2 = PBS_Fun(Fst_1, Fst_3, Fst_2)
#    PBS3 = PBS_Fun(Fst_3, Fst_2, Fst_1)
#    # median of Ts
#    Tmed = stats.median(PBS1[1:])
#    PBS_med = stats.median([PBS1[0],PBS2[0],PBS3[0]])
#    if Tmed > 0:
#        correction_term = PBS1[3]*PBS_med/Tmed
#    else:
#        correction_term = PBS1[3]
#    return(PBS1[0] - correction_term)

#  PBSN1 - another correction for long branches throughout the tree.
# Crawford et al 2016
def PBSN1_Fun(Fst_1, Fst_2, Fst_3):
    PBS1 = PBS_Fun(Fst_1, Fst_2, Fst_3)
    PBS2 = PBS_Fun(Fst_1, Fst_3, Fst_2)
    PBS3 = PBS_Fun(Fst_3, Fst_2, Fst_1)
    # median of Ts
    PBSN1 = PBS1[0]/(1 + PBS1[0] + PBS2[0] + PBS3[0])
    return(PBSN1)


samp_size = 25
# rescale for 50x
# specification of demographic history doesn't vary across replicates
demography = msprime.Demography()
demography.add_population(name = "A", initial_size = 4e+4, initially_active=True)
demography.add_population(name = "B", initial_size = 2e+4, initially_active=True)
demography.add_population(name = "C", initial_size = 2e+4)
demography.add_population_parameters_change(time=2e+3, population="B", initial_size=1e+4)
demography.add_population_split(time=2e+3, derived=["C"], ancestral="B")
demography.add_population_split(time=4e+3, derived=["B"], ancestral="A")


# rescale as 50x rather than 100x
# samp_size
def run_sim(num_samples, ancestry_seed, mutation_seed):
    ancestral_ts = msprime.sim_ancestry(samples={"A":samp_size, "B":samp_size, "C":samp_size}, demography = demography, sequence_length=5000, recombination_rate = 5.45e-7, gene_conversion_rate = 3.125e-6, gene_conversion_tract_length = 518, random_seed= ancestry_seed)    
    mutated_ts = msprime.sim_mutations(ancestral_ts, rate=2.605e-7, keep=True, random_seed=mutation_seed, model=msprime.SLiMMutationModel(type=0))  
    #mutated_ts = msprime.sim_mutations(ancestral_ts, rate=5.21e-8, keep=True, random_seed=mutation_seed, model=msprime.BinaryMutationModel())
    return mutated_ts

num_replicates = 1
# num_replicates = 10000
make_seed = random.randint(0,2**31)
rng = np.random.RandomState(make_seed)
seeds = rng.randint(1, 2**31, size=(num_replicates, 2))

outfile = open("ms_rep_output_dros", "w")

# no header row for single output files, since these are concatenated
# add header row afterwards

#outfile.write("Fst_AB" + '\t' + "Fst_AC" + '\t' + "Fst_BC" + '\t' + "PBS_A" + '\t' + "PBS_B" + '\t' + "PBS_C" + '\t' + "PBE_A" + '\t' + "PBE_B" + '\t' + "PBE_C" + '\t' + "PBSN1_A" + '\t' + "PBSN1_B" + '\t' + "PBSN1_C" + '\n')


#outfile.write("seg_sites_A" + '\t' + "seg_sites_B" + '\t' + "seg_sites_C" + '\t' + "pi_A" + '\t' + "pi_B" + '\t' + "pi_C" + '\t' + "Fst_AB" + '\t' + "Fst_AC" + '\t' + "Fst_BC" + '\t' + "PBS_A" + '\t' + "PBS_B" + '\t' + "PBS_C" + '\t' + "PBE_A" + '\t' + "PBE_B" + '\t' + "PBE_C" + '\t' + "PBSN1_A" + '\t' + "PBSN1_B" + '\t' + "PBSN1_C" + '\n')


#print(len(mts.samples(population=0)))
#print(len(mts.samples(population=1)))
#print(len(mts.samples(population=2)))

# treat C as focal population ?
for rep_index in range(num_replicates):
    mts = run_sim(num_replicates, *seeds[rep_index])
    seg_sites = mts.segregating_sites(sample_sets=[mts.samples(population=0), mts.samples(population=1), mts.samples(population=2),]).tolist()
    #seg_sites = seg_sites
    diversity = mts.diversity(sample_sets=[mts.samples(population=0), mts.samples(population=1), mts.samples(population=2),]).tolist()
    hap_size = 2*samp_size
    Genotype_Array = mts.genotype_matrix()
    Genotype_Array_A = Genotype_Array[:,:hap_size]
    Genotype_Array_B = Genotype_Array[:,hap_size:2*hap_size]
    Genotype_Array_C = Genotype_Array[:,2*hap_size:3*hap_size]
    site_freqs_a = np.sum(Genotype_Array_A, 1)/hap_size
    site_freqs_b = np.sum(Genotype_Array_B, 1)/hap_size
    site_freqs_c = np.sum(Genotype_Array_C, 1)/hap_size
    Fst_AB = weighted_avg_rey_fst(site_freqs_a.tolist(), site_freqs_b.tolist(), hap_size, hap_size)
    Fst_AC = weighted_avg_rey_fst(site_freqs_a.tolist(), site_freqs_c.tolist(), hap_size, hap_size)
    Fst_BC = weighted_avg_rey_fst(site_freqs_b.tolist(), site_freqs_c.tolist(), hap_size, hap_size)
    PBS_A = PBS_Fun(Fst_AB, Fst_AC, Fst_BC)[0]
    PBS_B = PBS_Fun(Fst_AB, Fst_BC, Fst_AC)[0]
    PBS_C = PBS_Fun(Fst_AC, Fst_BC, Fst_AB)[0]
    PBSN1_A = PBSN1_Fun(Fst_AB,Fst_AC, Fst_BC)
    PBSN1_B = PBSN1_Fun(Fst_AB,Fst_BC, Fst_AC)
    PBSN1_C = PBSN1_Fun(Fst_AC, Fst_BC, Fst_AB)
    outfile.write(str(seg_sites[0]) + '\t' + str(seg_sites[1]) + '\t' + str(seg_sites[2]) + '\t' + str(diversity[0]) + '\t' + str(diversity[1]) + '\t' + str(diversity[2]) + '\t' + str(Fst_AB) + '\t' + str(Fst_AC) + '\t' + str(Fst_BC) + '\t' + str(PBS_A) + '\t' + str(PBS_B) + '\t' + str(PBS_C) + '\t' + str(PBSN1_A) + '\t' + str(PBSN1_B) + '\t' + str(PBSN1_C) + '\n')
    
outfile.close()

#print sequence data
GA = np.transpose(Genotype_Array)

outfile = open("samp_seq_rep_output_dros", 'w')
for row in GA:
    outfile.write(' '.join([str(a) for a in row]) + '\n')    
