mkr<-read.table("/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_cells/z20191101_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s5_ctg_markers_stat.txt")

mkrall <- mkr[mkr$V3=="P&M-cluster",]

hist(mkrall$V7, breaks=1000, xlim=c(20, 150), xlab="marker coverage a long a contig", ylab="Number of contigs", main="Marker coverage", cex.main=0.8)

###
x895<-read.table("/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/reads_895cells")

y895<- x895[x895$V1<100000, ]

y895<- x895[x895$V1<200000, ]