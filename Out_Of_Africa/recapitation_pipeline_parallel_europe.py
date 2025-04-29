# recapitation of 4 population genealogy - Africa, out of Africa  with migration
# import Bvalue used in slim simulation
# Europe as focal population

import msprime, tskit, pyslim
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

#note on versions
#running pyslim 0.700, tskit 0.4.1, msprime 1.2.0

#ts = tskit.load("Japan_ThreePopFix.trees")
ts = tskit.load("OutOfAfrica_Parallel.trees")

# convert to format recognized by pyslim 0.700 (apparently !!!)
ts_new = ts
#ts_new = pyslim.SlimTreeSequence(ts)

num_replicates = 1
# num_replicates = 10000
make_seed = random.randint(0,2**31)
rng = np.random.RandomState(make_seed)
recap_seed = rng.randint(1,2**31)
ancestry_seed = rng.randint(1,2**31)
mutation_seed = rng.randint(1,2*31)

# recapitation
# rather than ancestral size, need to simulate demography - use only in 3 population case with no simulated splits

infile = open("Bval.txt", "r")
bval = infile.read().split(' ')[0].strip()
bval = float(bval)
#print(bval)

demography = msprime.Demography.from_tree_sequence(ts_new)
# tryp pop_1 for p1, what is p0??
demography.add_population_parameters_change(time=12310, population="p1", initial_size=int(7300*bval),growth_rate = 0.0)
#demography.add_population_parameters_change(time=12310, population="pop_1", initial_size=int(7300*bval),growth_rate = 0.0)

rts = pyslim.recapitate(ts_new, recombination_rate = 1.14e-8, gene_conversion_rate = 5.9e-8, gene_conversion_tract_length = 100, ancestral_Ne=int(7300*bval), random_seed = ancestry_seed)


# add mutations
#mts = msprime.mutate(rts, rate= 1.22e-8, keep=True, random_seed=10)
#mts = msprime.sim_mutations(rts, rate = 1.22e-8, keep=True, random_seed = mutation_seed, model=msprime.BinaryMutationModel())
mts = msprime.sim_mutations(rts, rate = 2.35e-8, keep=True, random_seed = mutation_seed, model=msprime.SLiMMutationModel(type=0))

outfile = open("slim_europe_rep_win_output", "w")
outfile2 = open("slim_europe_rep_site_output", "w")
outfile3 = open("slim_europe_rep_seq_a_output", "w")
outfile4 = open("slim_europe_rep_seq_b_output", "w")
outfile5 = open("slim_europe_rep_seq_c_output", "w")

samp_size = 25
hap_size = 2*samp_size
#mts = run_sim(num_replicates, *seeds[rep_index])
seg_sites = mts.segregating_sites(sample_sets=[mts.samples(population=1), mts.samples(population=2), mts.samples(population=3),]).tolist()
diversity = mts.diversity(sample_sets=[mts.samples(population=1), mts.samples(population=2), mts.samples(population=3),]).tolist()
mts_subset_a = mts.simplify(mts.samples(population=1), filter_sites=False)
mts_subset_b = mts.simplify(mts.samples(population=2), filter_sites=False)
mts_subset_c = mts.simplify(mts.samples(population=3), filter_sites=False)
Genotypes_A = mts_subset_a.genotype_matrix()
Genotypes_B = mts_subset_b.genotype_matrix()
Genotypes_C = mts_subset_c.genotype_matrix()
samp_a = np.random.choice(range(len(Genotypes_A[0])), size=hap_size, replace=False)
samp_b = np.random.choice(range(len(Genotypes_B[0])), size=hap_size, replace=False)
samp_c = np.random.choice(range(len(Genotypes_C[0])), size=hap_size, replace=False)
Genotype_Array_A = Genotypes_A[:,samp_a]
Genotype_Array_B = Genotypes_B[:,samp_b]
Genotype_Array_C = Genotypes_C[:,samp_c]
site_freqs_a = np.sum(Genotype_Array_A, 1)/hap_size
site_freqs_b = np.sum(Genotype_Array_B, 1)/hap_size
## for population C, check which site is near fixation
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
#GA = np.transpose(Genotype_Array)
#for row in GA:
#    outfile3.write(' '.join([str(a) for a in row]) + '\n')

GA = np.transpose(Genotype_Array_A)
for row in GA:
    outfile3.write(' '.join([str(a) for a in row]) + '\n')
GB = np.transpose(Genotype_Array_B)
for row in GB:
    outfile4.write(' '.join([str(a) for a in row]) + '\n')
GC = np.transpose(Genotype_Array_C)
for row in GC:
    outfile5.write(' '.join([str(a) for a in row]) + '\n')

outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()
outfile5.close()
