# simulate out of Africa demogrpahic model from Gutenkunst et al 2009

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
    
# for site fst, need to correct for 0,1 case with pseudocount = 0.5
def site_rey_fst(freq1, freq2, SampleSize1, SampleSize2):
    n1 = SampleSize1
    n2 = SampleSize2
    if freq1 < 0:
        freq1 = 0
    if freq2 < 0:
        freq2 = 0
    if freq1 > 1:
        freq1 = 1
    if freq2 > 1:
        freq2 = 1
    if (freq1 == 0 and freq2 == 1):
        freq1 = 1/(2*n1)
    if (freq1 == 1 and freq2 == 0):
        freq2 = 1/(2*n2)
    #print(freq1)
    #print(freq2)
    #print('\n')
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
    if (FST > 0.979):
        FST = 0.9791666666666667
    #if (FST < 0):
       # FST = 0
    #else:
        #FST = FST
    return([FST, WholeNum, WholeDen])
    
#

# calculate numerator and denominator terms of Reynolds Fst for every locus          
def multilocus_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
    fst_list = [rey_fst(Freqlist1[i], Freqlist2[i], SampleSize1, SampleSize2) for i in range(len(Freqlist1))]
    return(fst_list)


#maximum site fst in a window
def max_fst(Freqlist1, Freqlist2, SampleSize1, SampleSize2):
    fst_list = [site_rey_fst(Freqlist1[i], Freqlist2[i], SampleSize1, SampleSize2)[0] for i in range(len(Freqlist1))]
    maximum_fst = max(fst_list)
    return(maximum_fst)
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
        temp_ratio = 1.0
    else:
        temp_ratio = temp_ratio
    return(temp_ratio)
    
# PBS = Population Branch Statistic (Yi et al 2010)
# Fst1: AB, Fst2: AC, Fst3: BC
# T: rescale as -log(1-Fst) for approximate additivity
# PBS branch length estimate, using same methodology as NJ   
def PBS_Fun(Fst_1, Fst_2, Fst_3):
    #print(Fst_1)
    T1 = -math.log(1-Fst_1)
    T2 = -math.log(1-Fst_2)
    T3 = -math.log(1-Fst_3)
    PBS = (T1 + T2 - T3)/2
    return([PBS,T1,T2,T3])
    
# maximum site PBS in window
#def PBS_Max(Freqlist_1, Freqlist_2, Freqlist_3, SampleSize1, SampleSize2, SampleSize3):
#   fst_list_12 = [rey_fst(Freqlist1[i], Freqlist2[i], SampleSize1, SampleSize2) for i in range(len(Freqlist1))]
#    fst_list_13 = [rey_fst(Freqlist1[i], Freqlist3[i], SampleSize1, SampleSize3) for i in range(len(Freqlist1))]
#    fst_list_23 = [rey_fst(Freqlist2[i], Freqlist3[i], SampleSize2, SampleSize3) for i in range(len(Freqlist2))]
#    T1 = - math.log(1-fst_list_12)
#    T2 = - math.log(1-fst_list_13)
#    T3 = - math.log(1-fst_list_23)
#    PBS = (T1 + T2 - T3)/2
#    max_pbs = max(PBS)
#    return(max_pbs)

#def Replace_If_Neg(x):
#    if x <= 0:
#        return(0.0000001)
#    else:
#       return(x)
    
# PBS site max: PBS_B and PBS_C at site maximum of PBS_A as well as max T_BC outgroup
# use to calculate PBS site max (with maximum PBS per window fixed by site), as well as max T_BC for PBE calculation
def PBS_Site_Max(Freqlist_1, Freqlist_2, Freqlist_3, SampleSize1, SampleSize2, SampleSize3, MedWinPBS1, MedWinT23):
    fst_list_12 = [site_rey_fst(Freqlist_1[i], Freqlist_2[i], SampleSize1, SampleSize2)[0] for i in range(len(Freqlist_1))]
    fst_list_13 = [site_rey_fst(Freqlist_1[i], Freqlist_3[i], SampleSize1, SampleSize3)[0] for i in range(len(Freqlist_1))]
    fst_list_23 = [site_rey_fst(Freqlist_2[i], Freqlist_3[i], SampleSize2, SampleSize3)[0] for i in range(len(Freqlist_2))]
    #print(fst_list_12)
    fst_12_diff = [1-x for x in fst_list_12]
    fst_13_diff = [1-x for x in fst_list_13]
    fst_23_diff = [1-x for x in fst_list_23]
    #fst_12_diff = list(map(Replace_If_Neg, fst_12_diff))
    #fst_13_diff = list(map(Replace_If_Neg, fst_13_diff))
    #fst_23_diff = list(map(Replace_If_Neg, fst_13_diff))
    #print(min(fst_12_diff))
    #print(min(fst_13_diff))
    #print(min(fst_23_diff))
    T12 = list(map(math.log, fst_12_diff))
    T12 = [-x for x in T12]
    T13 = list(map(math.log, fst_13_diff))
    T13 = [-x for x in T13]
    T23 = list(map(math.log, fst_23_diff))
    T23 = [-x for x in T23]
    PBS_1 = [(T12[i] + T13[i] - T23[i])/2 for i in range(len(T12))]
    PBS_2 = [(T12[i] + T23[i] - T13[i])/2 for i in range(len(T12))]
    PBS_3 = [(T23[i] + T13[i] - T12[i])/2 for i in range(len(T23))]
    max_pbs = max(PBS_1)
    posits = [i for i in range(len(PBS_1)) if PBS_1[i] == max_pbs]
    #print(len(posits))
    pos = posits[0]
    # this is used to find PBS - T maximal for use in PBE*
    PBS_Diff_T1 = [PBS_1[i] - T23[i]*MedWinPBS1/MedWinT23 for i in range(len(T23))]
    max_pbs_tdiff = max(PBS_Diff_T1)
    posits2 = [i for i in range(len(PBS_Diff_T1)) if PBS_Diff_T1[i] == max_pbs_tdiff]
    pos2 = posits2[0]
    # note - T23 is outgroup to 1, so  it is in the same output as max(PBS_1) etc
    return([max_pbs, PBS_2[pos], PBS_3[pos], T23[pos], PBS_1[pos2], T23[pos2]])
    
def PBSN1_SiteMax(Freqlist_1, Freqlist_2, Freqlist_3, SampleSize1, SampleSize2, SampleSize3):
    dat_vec = PBS_Site_Max(Freqlist_1, Freqlist_2, Freqlist_3, SampleSize1, SampleSize2, SampleSize3)
    PBSN1_Max = dat_vec[0]/(1 + dat_vec[0] + dat_vec[1] + dat_vec[2])
    return(PBSN1_Max)
    
    
# maximum rescaled branch lengt T = - log(1-fst)
def T_Max(Freqlist_1, Freqlist_2, SampleSize1, SampleSize2):
    Fst = [site_rey_fst(Freqlist_1[i], Freqlist_2[i], SampleSize1, SampleSize2)[0] for i in range(len(Freqlist_1))]
    Fst_Diff = [1-x for x in Fst]
    #Fst_Diff = list(map(Replace_If_Neg, Fst_Diff))
    T1 = list(map(math.log, Fst_Diff))
    T1 = [-x for x in T1]
    max_t = max(T1)
    return(max_t)

#  PBSN1 - another correction for long branches throughout the tree.
# Crawford et al 2016
def PBSN1_Fun(Fst_1, Fst_2, Fst_3):
    PBS1 = PBS_Fun(Fst_1, Fst_2, Fst_3)
    PBS2 = PBS_Fun(Fst_1, Fst_3, Fst_2)
    PBS3 = PBS_Fun(Fst_3, Fst_2, Fst_1)
    # median of Ts
    PBSN1 = PBS1[0]/(1 + PBS1[0] + PBS2[0] + PBS3[0])
    return(PBSN1)


# obtained from previous run without PBE*, used as input
PBS_A_Win_Med = 0.1221635
PBS_B_Win_Med = 0.0001883063
PBS_C_Win_Med = 0.06006449
T_AB_Win_Med = 0.1270
T_AC_Win_Med = 0.1903012
T_BC_Win_Med = 0.06256379

Bvals = [0.4145932, 0.4447721, 0.4747853, 0.5050735, 0.5343682, 0.5646864, 0.5948979, 0.6247126, 0.6544159, 0.6849211, 0.7148560, 0.7450127, 0.7751756, 0.8053146, 0.8351195, 0.8650366, 0.8949025, 0.9247594, 0.9546144, 0.9872088]

Bval_Freq = [0.01239254, 0.02456442, 0.03853131, 0.05426460, 0.07071475, 0.08913901, 0.11009095, 0.13321839, 0.16011891, 0.18961042, 0.22439077, 0.26000370, 0.30348732, 0.35492159, 0.41352148, 0.48974020, 0.58184641, 0.69745901, 0.83316061, 1.00000000]

# Gutenkunst model mutation rate 2.358e-8, vs. 1.14e-8 in our paper, 1.29e-8 in other publications
def run_sim(num_samples, ancestry_seed, mutation_seed):
    ancestral_ts = msprime.sim_ancestry(samples={"A":samp_size, "B":samp_size, "C":samp_size}, demography= demography, sequence_length=100000, recombination_rate=1.14e-8, gene_conversion_rate=5.9e-8, gene_conversion_tract_length=100, random_seed= ancestry_seed)    
    mutated_ts = msprime.sim_mutations(ancestral_ts, rate=2.35e-8, keep=True, random_seed=mutation_seed, model=msprime.SLiMMutationModel(type=0))  
    #mutated_ts = msprime.sim_mutations(ancestral_ts, rate=5.21e-8, keep=True, random_seed=mutation_seed, model=msprime.BinaryMutationModel())
    return mutated_ts
    
num_replicates = 1000000
#num_replicates = 100
make_seed = random.randint(0,2**31)
rng = np.random.RandomState(make_seed)
seeds = rng.randint(1, 2**31, size=(num_replicates, 2))

outfile = open("new_demog_window_europe", "w")
outfile2 = open("new_demog_site_ms_europe", "w")
outfile3 = open("new_demog_samp_seq_europe", "w")

# treat C as focal population (FST_BC, PBS_C etc), T_AB as outgroup branch for PBE
for rep_index in range(num_replicates):

    pick_bin = np.random.uniform()
    compare_probs = [el for el in Bval_Freq if el <= pick_bin]
    pos = len(compare_probs)
    sim_b = Bvals[pos]

    samp_size = 25
    # specification of demographic history doesn't vary across replicates
    demography = msprime.Demography()
    demography.add_population(name = "A", initial_size = 12300*sim_b, initially_active=True)
    demography.add_population(name = "B", initial_size = 295725*sim_b, initially_active=True, growth_rate = 0.004)
    demography.add_population(name = "C", initial_size = 53404*sim_b, growth_rate = 0.0055)

    # model 0 migration by commenting out set_migration_rate and add_migration_rate_change
    #demography.set_symmetric_migration_rate(["A","B","C"],0.0001)
    demography.set_symmetric_migration_rate(["A","B"], rate = 0.00003)
    demography.set_symmetric_migration_rate(["A","C"], rate = 0.000019)
    demography.set_symmetric_migration_rate(["B","C"], rate = 0.000096)
    demography.add_migration_rate_change(time = 848, rate = 0.00025, source = "B", dest = "A")
    demography.add_migration_rate_change(time = 848, rate = 0.00025, source = "A", dest = "B")

    demography.add_population_parameters_change(time=848, population="B", initial_size=2100*sim_b,growth_rate = 0.0)
    #demography.add_population_parameters_change(time=8800, population="A", initial_size=7300*sim_b, growth_rate = 0.0)
    demography.add_population_split(time=848, derived=["C"], ancestral="B")
    demography.add_population_split(time=5600, derived=["B"], ancestral="A")
    demography.add_population_parameters_change(time=8800, population="A", initial_size=7300*sim_b, growth_rate = 0.0)

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
    
    Fst_SiteMax_AB = max_fst(site_freqs_a.tolist(), site_freqs_b.tolist(), hap_size, hap_size)
    Fst_SiteMax_AC = max_fst(site_freqs_a.tolist(), site_freqs_c.tolist(), hap_size, hap_size)
    Fst_SiteMax_BC = max_fst(site_freqs_b.tolist(), site_freqs_c.tolist(), hap_size, hap_size)
    
    #PBS_Site_Max(Freqlist_1, Freqlist_2, Freqlist_3, SampleSize1, SampleSize2, SampleSize3, MedWinPBS1, MedWinT23):
    PBS_SiteMax_Stats_A = PBS_Site_Max(site_freqs_a.tolist(), site_freqs_b.tolist(), site_freqs_c.tolist(), hap_size, hap_size, hap_size, PBS_A_Win_Med, T_BC_Win_Med)
    PBS_SiteMax_Stats_B = PBS_Site_Max(site_freqs_b.tolist(), site_freqs_a.tolist(), site_freqs_c.tolist(), hap_size, hap_size, hap_size, PBS_B_Win_Med, T_AC_Win_Med)
    PBS_SiteMax_Stats_C = PBS_Site_Max(site_freqs_c.tolist(), site_freqs_a.tolist(), site_freqs_b.tolist(), hap_size, hap_size, hap_size, PBS_C_Win_Med, T_AB_Win_Med)
    
    PBS_SiteMax_A = PBS_SiteMax_Stats_A[0]
    PBS_SiteMax_B = PBS_SiteMax_Stats_B[0]
    PBS_SiteMax_C = PBS_SiteMax_Stats_C[0]
    
    PBSN1_SiteMax_A = PBS_SiteMax_Stats_A[0]/(1 + PBS_SiteMax_Stats_A[0] + PBS_SiteMax_Stats_A[1] + PBS_SiteMax_Stats_A[2])
    PBSN1_SiteMax_B = PBS_SiteMax_Stats_B[0]/(1 + PBS_SiteMax_Stats_B[0] + PBS_SiteMax_Stats_B[1] + PBS_SiteMax_Stats_B[2])
    PBSN1_SiteMax_C = PBS_SiteMax_Stats_C[0]/(1 + PBS_SiteMax_Stats_C[0] + PBS_SiteMax_Stats_C[1] + PBS_SiteMax_Stats_C[2])
     
    # for use in PBE* calculations, T_AB is outgroup to C etc
    TStar_AB = PBS_SiteMax_Stats_C[5]
    TStar_AC = PBS_SiteMax_Stats_B[5]
    TStar_BC = PBS_SiteMax_Stats_A[5]
    
    # PBS for use in PBE* calculations
    PBS_PBE_Star_A = PBS_SiteMax_Stats_A[4]
    PBS_PBE_Star_B = PBS_SiteMax_Stats_B[4]
    PBS_PBE_Star_C = PBS_SiteMax_Stats_C[4]
    
    # for use in PBE calculations
    TMax_AB = T_Max(site_freqs_a.tolist(), site_freqs_b.tolist(), hap_size, hap_size)
    TMax_AC = T_Max(site_freqs_a.tolist(), site_freqs_c.tolist(), hap_size, hap_size)
    TMax_BC = T_Max(site_freqs_b.tolist(), site_freqs_c.tolist(), hap_size, hap_size)
    
    outfile.write(str(seg_sites[0]) + '\t' + str(seg_sites[1]) + '\t' + str(seg_sites[2]) + '\t' + str(diversity[0]) + '\t' + str(diversity[1]) + '\t' + str(diversity[2]) + '\t' + str(Fst_AB) + '\t' + str(Fst_AC) + '\t' + str(Fst_BC) + '\t' + str(PBS_A) + '\t' + str(PBS_B) + '\t' + str(PBS_C) + '\t' + str(PBSN1_A) + '\t' + str(PBSN1_B) + '\t' + str(PBSN1_C) + '\n')
    
    outfile2.write(str(Fst_SiteMax_AB) + '\t' + str(Fst_SiteMax_AC) + '\t' + str(Fst_SiteMax_BC) + '\t' + str(PBS_SiteMax_A) + '\t' + str(PBS_SiteMax_B) + '\t' + str(PBS_SiteMax_C) + '\t' + str(PBSN1_SiteMax_A) + '\t' + str(PBSN1_SiteMax_B) + '\t' + str(PBSN1_SiteMax_C) + '\t' + str(TStar_AB) + '\t' + str(TStar_AC) + '\t' + str(TStar_BC) + '\t' + str(PBS_PBE_Star_A) + '\t' + str(PBS_PBE_Star_B) + '\t' + str(PBS_PBE_Star_C) + '\t' + str(TMax_AB) + '\t' + str(TMax_AC) + '\t' + str(TMax_BC) + '\n')
    
    GA = np.transpose(Genotype_Array)
    for row in GA:
        outfile3.write(' '.join([str(a) for a in row]) + '\n')
    #outfile.write(str(Fst_AB) + '\t' + str(Fst_AC) + '\t' + str(Fst_BC) + '\t' + str(PBS_A) + '\t' + str(PBS_B) + '\t' + str(PBS_C) + '\t' + str(PBE_A) + '\t' + str(PBE_B) + '\t' + str(PBE_C) + '\t' + str(PBSN1_A) + '\t' + str(PBSN1_B) + '\t' + str(PBSN1_C) + '\n')

outfile.close()
outfile2.close()
outfile3.close()
