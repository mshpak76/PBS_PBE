# for the case where the full 2K generation genealogy is simulated in slim
# based on SlimCode_Fixation_Revised_ToRecap.txt code
# recapitate "prehistory" - no need for input demography, only Ne

# recapitates tree from three population model (1000 generations) and no splits
# calculates Fst, PBS, PBSn1 from output

# includes migration such that Nm=1 from each population


import msprime, tskit, pyslim
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

#note on versions
#running pyslim 0.700, tskit 0.4.1, msprime 1.2.0

#ts = tskit.load("OnePop.trees")
ts = tskit.load("ThreePopFixMigration.trees")

# convert to format recognized by pyslim 0.700 (apparently !!!)
ts_new = pyslim.SlimTreeSequence(ts)

num_replicates = 1
# num_replicates = 10000
make_seed = random.randint(0,2**31)
rng = np.random.RandomState(make_seed)
recap_seed = rng.randint(1,2**31)
ancestry_seed = rng.randint(1,2**31)
mutation_seed = rng.randint(1,2*31)

# recapitation
# rather than ancestral size, need to simulate demography - use only in 3 population case with no simulated splits

rts = pyslim.recapitate(ts_new, recombination_rate = 1.14e-8, gene_conversion_rate = 5.9e-8, gene_conversion_tract_length = 100, ancestral_Ne=20000, random_seed = ancestry_seed)


# add mutations
#mts = msprime.mutate(rts, rate= 1.22e-8, keep=True, random_seed=10)

#mts = msprime.sim_mutations(rts, rate = 1.22e-8, keep=True, random_seed = mutation_seed, model=msprime.BinaryMutationModel())

mts = msprime.sim_mutations(rts, rate = 1.22e-8, keep=True, random_seed = mutation_seed, model=msprime.SLiMMutationModel(type=0))

outfile = open("slim_sweep_migration_output", "w")
outfile2 = open("sweep_mig_seq_a", "w")
outfile3 = open("sweep_mig_seq_b", "w")
outfile4 = open("sweep_mig_seq_c", "w")

samp_size = 25
hap_size = samp_size * 2
seg_sites = mts.segregating_sites(sample_sets=[mts.samples(population=1), mts.samples(population=2), mts.samples(population=3),]).tolist()
diversity = mts.diversity(sample_sets=[mts.samples(population=1), mts.samples(population=2), mts.samples(population=3),]).tolist()
mts_subset_a = mts.simplify(mts.samples(population=1), filter_sites=False)
mts_subset_b = mts.simplify(mts.samples(population=2), filter_sites=False)
mts_subset_c = mts.simplify(mts.samples(population=3), filter_sites=False)
Genotype_Array = mts.genotype_matrix()
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
#PBE_A = PBE_Fun(Fst_AB,Fst_AC, Fst_BC)
#PBE_B = PBE_Fun(Fst_AB,Fst_BC, Fst_AC)
#PBE_C = PBE_Fun(Fst_AC, Fst_BC, Fst_AB)
PBSN1_A = PBSN1_Fun(Fst_AB,Fst_AC, Fst_BC)
PBSN1_B = PBSN1_Fun(Fst_AB,Fst_BC, Fst_AC)
PBSN1_C = PBSN1_Fun(Fst_AC, Fst_BC, Fst_AB)
#outfile.write("seg_sites_A" + "\t" + "seg_sites_B" + "\t" + "seg_sites_C" + "\t" + "nuc_diversity_A" + "\t" + "nuc_diversity_B" + "\t" + "nuc_diversity_C" + "\t" + "Fst_AB" + "\t" + "Fst_AC" + "\t" + "Fst_BC" + "\t" + "PBS_A" + "\t" + "PBS_B" + "\t" + "PBS_C" + "\t" +  "PBSN1_A" + "\t" + "PBSN1_B" + "\t" + "PBSN1_C" + '\n')
outfile.write(str(seg_sites[0]) + '\t' + str(seg_sites[1]) + '\t' + str(seg_sites[2]) + '\t' + str(diversity[0]) + '\t' + str(diversity[1]) + '\t' + str(diversity[2]) + '\t' + str(Fst_AB) + '\t' + str(Fst_AC) + '\t' + str(Fst_BC) + '\t' + str(PBS_A) + '\t' + str(PBS_B) + '\t' + str(PBS_C) + '\t' + str(PBSN1_A) + '\t' + str(PBSN1_B) + '\t' + str(PBSN1_C) + '\n')
GA = np.transpose(Genotype_Array_A)
for row in GA:
    outfile2.write(' '.join([str(a) for a in row]) + '\n')
GB = np.transpose(Genotype_Array_B)
for row in GB:
    outfile3.write(' '.join([str(a) for a in row]) + '\n')
GC = np.transpose(Genotype_Array_C)
for row in GC:
    outfile4.write(' '.join([str(a) for a in row]) + '\n')
    #outfile.write(str(Fst_AB) + '\t' + str(Fst_AC) + '\t' + str(Fst_BC) + '\t' + str(PBS_A) + '\t' + str(PBS_B) + '\t' + str(PBS_C) + '\t' + str(PBE_A) + '\t' + str(PBE_B) + '\t' + str(PBE_C) + '\t' + str(PBSN1_A) + '\t' + str(PBSN1_B) + '\t' + str(PBSN1_C) + '\n')

outfile.close()
outfile2.close()
outfile3.close()
outfile4.close()
