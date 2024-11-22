# create simulated "genomes"
# 98% neutral (or bgs), 1% local, 2% shared
# 99% neutral, 1% local

setwd("/raid10/mshpak/Fst_Simulations/SlimScripts/Run_MSPrime/Run_MSPrime_CHTC/Sequence_Data/Neutral_Dros")
neutral = read.csv("site_msprime_output_largepop.csv", header=TRUE, sep='\t')

setwd("/raid10/mshpak/Fst_Simulations/SlimScripts/Recapitation/Soft/Sequences/Dros_005/")
local_sweep = read.csv("slim_site_fst_output.csv",header=TRUE, sep='\t')
setwd("/raid10/mshpak/Fst_Simulations/SlimScripts/Recapitation/Parallel/Dros/Sequences/Soft_005/")
global_sweep = read.csv("slim_site_fst_output.csv",header=TRUE,sep='\t')

# need window pbs, fst for PBE calculations
setwd("/raid10/mshpak/Fst_Simulations/SlimScripts/Recapitation/Analyze_Genomes")
neutral_win = read.csv("merged_neutral_largepop.csv", header=TRUE, sep='\t')
local_sweep_win = read.csv("sweep_soft005_largepop_05.csv", header=TRUE, sep='\t')
global_sweep_win  = read.csv("merge_sweep05_parallel_soft005_largepop.csv", header=TRUE, sep='\t')

neutral = signif(neutral,7)
local_sweep = signif(local_sweep,7)
global_sweep = signif(global_sweep,7)

#local_sweep = read.csv("slim_sweep05_part80_largepop.csv", header=TRUE, sep='\t')
#local_sweep = read.csv("slim_sweep05_part50_largepop.csv", header=TRUE, sep='\t')
#background = read.csv("merged_bgs_largepop.csv", header=TRUE, sep='\t')
#global_sweep = read.csv("slim_threesweep_smallpop_025.csv", header=TRUE, sep='\t')
#global_sweep  = read.csv("parallel05_partial50_largepop_merged.csv", header=TRUE, sep='\t')
#global_sweep = read.csv("slim_sweep025_parallel_part80_smallpop.csv", header=TRUE, sep='\t')
#global_sweep = read.csv("slim_parallel05_partial80_largepop.csv", header=TRUE, sep='\t')


nloci = 25000
# fraction neutral vs. local selection, global selection
#sel_freq = (0.98, 0.01, 0.01)
# alternative scenario 99% neutral
sel_freq = c(0.99,0.01)


# substitute bgs for local in both versions
# neutral or local selection
nreps = 1000 #number of replicates with model genome to test, start with tractable number

# statistics: Fst_BC, PBS_C, PBE_C, PBSN1_C
# calculate PBE from PBS etc

setwd("/raid10/mshpak/Fst_Simulations/SlimScripts/Recapitation/Analyze_Genomes/Site_Fst")

# true positive, false positive, true negative, false negative
fst = c()
pbs = c()
pbsn1 = c()
pbe = c()
quant_pval = 0.99
for(i in 1:nreps){
	nsel = rbinom(1,nloci,0.01)
	nneut = nloci - nsel
	sample_neut = sample(1:1000000, nneut, replace = FALSE)
	sample_sel = sample(1:10000, nsel, replace = FALSE)
	model_genome = rbind(neutral[sample_neut,], local_sweep[sample_sel,])
	# for PBE - need window statistics
	model_genome_win = rbind(neutral_win[sample_neut,],local_sweep_win[sample_sel,])
	# to calculate PBE
	PBS_Med_C = median(model_genome_win$PBS_C)
	T_AB = -log(1-model_genome_win$Fst_AB)
	PBS_Med_B = median(model_genome_win$PBS_B)
	T_AC = -log(1-model_genome_win$Fst_AC)
	PBS_Med_A = median(model_genome_win$PBS_A)
	T_BC = -log(1-model_genome_win$Fst_BC)
	# use window PBS and T medians, compare to site PBS
	PBE_C = model_genome$PBS_C - T_AB*PBS_Med_C/median(T_AB)
	PBE_B = model_genome$PBS_B - T_AC*PBS_Med_B/median(T_AC)
	PBE_A = model_genome$PBS_A - T_BC*PBS_Med_A/median(T_BC)
	newmat = cbind(cbind(cbind(model_genome, PBE_C),PBE_B),PBE_A)
	#99% quantiles
	model_genome = newmat
	fstbc_quant01 = quantile(model_genome$Fst_BC, quant_pval)
	pbsc_quant01 = quantile(model_genome$PBS_C, quant_pval)
	pbsn1c_quant01 = quantile(model_genome$PBSN1_C, quant_pval)
	pbec_quant01 = quantile(model_genome$PBE_C, quant_pval)
	fst_top = which(newmat$Fst_BC >= fstbc_quant01)
	pbs_top = which(newmat$PBS_C >= pbsc_quant01)
	pbsn1_top = which(newmat$PBSN1_C >= pbsn1c_quant01)
	pbe_top = which(newmat$PBE_C >= pbec_quant01)
	#which are true vs. false positive according to quantile criterion
	#pos1 = round(nloci*(1-quant_pval))
	pos_fst = length(fst_top)
	neg_fst = nloci - pos_fst
	pos_pbs = length(pbs_top)
	neg_pbs = nloci - pos_pbs
	pos_pbsn1 = length(pbsn1_top)
	neg_pbsn1 = nloci - pos_pbsn1
	pos_pbe = length(pbe_top)
	neg_pbe = nloci - pos_pbe
	
	fst_true_pos = length(which(fst_top > nneut))
	fst_false_pos = pos_fst - fst_true_pos 
	fst_false_neg = nsel - fst_true_pos
	fst_true_neg = neg_fst - fst_false_neg
	fst = rbind(fst, c(fst_true_pos, fst_false_pos, fst_true_neg, fst_false_neg ))
	colnames(fst) = c("true_pos", "false_pos", "true_neg", "false_neg")
	
	pbs_true_pos = length(which(pbs_top > nneut))
	pbs_false_pos = pos_pbs - pbs_true_pos
	pbs_false_neg = nsel - pbs_true_pos
	pbs_true_neg = neg_pbs - pbs_false_neg
	pbs = rbind(pbs, c(pbs_true_pos, pbs_false_pos, pbs_true_neg, pbs_false_neg))
	colnames(pbs) = c("true_pos", "false_pos", "true_neg", "false_neg")
	
	pbsn1_true_pos = length(which(pbsn1_top > nneut))
	pbsn1_false_pos = pos_pbsn1 - pbsn1_true_pos
	pbsn1_false_neg = nsel - pbsn1_true_pos
	pbsn1_true_neg = neg_pbsn1 - pbsn1_false_neg
	pbsn1 = rbind(pbsn1, c(pbsn1_true_pos, pbsn1_false_pos, pbsn1_true_neg, pbsn1_false_neg))
	colnames(pbsn1) = c("true_pos", "false_pos", "true_neg", "false_neg")
	
	pbe_true_pos = length(which(pbe_top > nneut))
	pbe_false_pos = pos_pbe - pbe_true_pos
	pbe_false_neg = nsel - pbe_true_pos
	pbe_true_neg = neg_pbe - pbe_false_neg
	pbe = rbind(pbe, c(pbe_true_pos, pbe_false_pos, pbe_true_neg, pbe_false_neg))
	colnames(pbe) = c("true_pos", "false_pos", "true_neg", "false_neg")
	#summary statistics comparison here
}

#

write.table(fst, file="Neutral_Sweep_Local_Soft005_Largepop_Site_Fst.csv", row.names=FALSE, col.names=TRUE,sep='\t')
write.table(pbs, file="Neutral_Sweep_Local_Soft005_Largepop_Site_PBS.csv", row.names=FALSE, col.names=TRUE,sep='\t')
write.table(pbsn1, file="Neutral_Sweep_Local_Soft005_Largepop_PBSN1.csv", row.names=FALSE, col.names=TRUE,sep='\t')
write.table(pbe, file="Neutral_Sweep_Local_Soft005_Largepop_PBE.csv", row.names=FALSE, col.names=TRUE,sep='\t')

# neutral or local or global
fst_2 = c()
pbs_2 = c()
pbsn1_2 = c()
pbe_2 = c()
quant_pval = 0.99
for(i in 1:nreps){
	class_counts = rmultinom(1, 25000,c(0.98,0.01,0.01))[,1]
	nneut = class_counts[1]
	global_sel = class_counts[2]
	local_sel = class_counts[3]
	not_local = nneut + global_sel
	sample_neut = sample(1:1000000, nneut, replace = FALSE)
	sample_global = sample(1:length(global_sweep[,1]), global_sel, replace = FALSE)
	sample_local = sample(1:10000, local_sel, replace = FALSE)
	model_genome = rbind(neutral[sample_neut,], global_sweep[sample_global,])
	model_genome = rbind(model_genome, local_sweep[sample_local,])
	# for pbe
	model_genome_win = rbind(neutral_win[sample_neut,], global_sweep_win[sample_global,])
	model_genome_win = rbind(model_genome_win, local_sweep_win[sample_local,])
	#calculate PBE
	PBS_Med_C = median(model_genome_win$PBS_C)
	T_AB = -log(1-model_genome_win$Fst_AB)
	PBS_Med_B = median(model_genome_win$PBS_B)
	T_AC = -log(1-model_genome_win$Fst_AC)
	PBS_Med_A = median(model_genome_win$PBS_A)
	T_BC = -log(1-model_genome_win$Fst_BC)
	PBE_C = model_genome$PBS_C - T_AB*PBS_Med_C/median(T_AB)
	PBE_B = model_genome$PBS_B - T_AC*PBS_Med_B/median(T_AC)
	PBE_A = model_genome$PBS_A - T_BC*PBS_Med_A/median(T_BC)
	newmat = cbind(cbind(cbind(model_genome, PBE_C),PBE_B),PBE_A)
	#99% quantiles
	model_genome = newmat
	fstbc_quant01 = quantile(model_genome$Fst_BC, quant_pval)
	pbsc_quant01 = quantile(model_genome$PBS_C, quant_pval)
	pbsn1c_quant01 = quantile(model_genome$PBSN1_C, quant_pval)
	pbec_quant01 = quantile(model_genome$PBE_C, quant_pval)
	fst_top = which(newmat$Fst_BC >= fstbc_quant01)
	pbs_top = which(newmat$PBS_C >= pbsc_quant01)
	pbsn1_top = which(newmat$PBSN1_C >= pbsn1c_quant01)
	pbe_top = which(newmat$PBE_C >= pbec_quant01)
	#which are true vs. false positive according to quantile criterion
	#pos1 = round(nloci*(1-quant_pval))
	pos_fst = length(fst_top)
	neg_fst = nloci - pos_fst
	pos_pbs = length(pbs_top)
	neg_pbs = nloci - pos_pbs
	pos_pbsn1 = length(pbsn1_top)
	neg_pbsn1 = nloci - pos_pbsn1
	pos_pbe = length(pbe_top)
	neg_pbe = nloci - pos_pbe
	
	fst_true_pos = length(which(fst_top > not_local))
	fst_false_pos = pos_fst - fst_true_pos 
	fst_false_neg = local_sel - fst_true_pos
	fst_true_neg = neg_fst - fst_false_neg
	fst_2 = rbind(fst_2, c(fst_true_pos, fst_false_pos, fst_true_neg, fst_false_neg ))
	colnames(fst_2) = c("true_pos", "false_pos", "true_neg", "false_neg")
	
	pbs_true_pos = length(which(pbs_top > not_local))
	pbs_false_pos = pos_pbs - pbs_true_pos
	pbs_false_neg = local_sel - pbs_true_pos
	pbs_true_neg = neg_pbs - pbs_false_neg
	pbs_2 = rbind(pbs_2, c(pbs_true_pos, pbs_false_pos, pbs_true_neg, pbs_false_neg))
	colnames(pbs_2) = c("true_pos", "false_pos", "true_neg", "false_neg")
	
	pbsn1_true_pos = length(which(pbsn1_top > not_local))
	pbsn1_false_pos = pos_pbsn1 - pbsn1_true_pos
	pbsn1_false_neg = local_sel - pbsn1_true_pos
	pbsn1_true_neg = neg_pbsn1 - pbsn1_false_neg
	pbsn1_2 = rbind(pbsn1_2, c(pbsn1_true_pos, pbsn1_false_pos, pbsn1_true_neg, pbsn1_false_neg))
	colnames(pbsn1_2) = c("true_pos", "false_pos", "true_neg", "false_neg")
	
	pbe_true_pos = length(which(pbe_top > not_local))
	pbe_false_pos = pos_pbe - pbe_true_pos
	pbe_false_neg = local_sel - pbe_true_pos
	pbe_true_neg = neg_pbe - pbe_false_neg
	pbe_2 = rbind(pbe_2, c(pbe_true_pos, pbe_false_pos, pbe_true_neg, pbe_false_neg))
	colnames(pbe_2) = c("true_pos", "false_pos", "true_neg", "false_neg")
	#summary statistics comparison here
}
#
write.table(fst_2, file="Neutral_Sweep_Global_Soft005_Largepop_Site_Fst.csv", row.names=FALSE, col.names=TRUE,sep='\t')
write.table(pbs_2, file="Neutral_Sweep_Global_Soft005_Largepop_Site_PBS.csv", row.names=FALSE, col.names=TRUE,sep='\t')
write.table(pbsn1_2, file="Neutral_Sweep_Global_Soft005_Largepop_Site_PBSN1.csv", row.names=FALSE, col.names=TRUE,sep='\t')
write.table(pbe_2, file="Neutral_Sweep_Global_Soft005_Largepop_Site_PBE.csv", row.names=FALSE, col.names=TRUE,sep='\t')
