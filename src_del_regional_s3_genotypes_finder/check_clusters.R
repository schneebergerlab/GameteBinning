x<-read.table("/biodata/dep_coupland/grp_schneeberger/projects/methods/src_shq/z_asPollinator_2019/src_genotype_del_regional_markers/test_tmp_pollen_del_like_genotypes_addi/s0_genotype_pollen_seq_del_wise/tmptorm.txt")

y<-t(x[x$V1=="tig00000013_pilon", ])

z<-y[4:length(y[, 1]), ]

k<- rowSums(z>0)

fit <- kmeans(k, 2) # 5 cluster solution