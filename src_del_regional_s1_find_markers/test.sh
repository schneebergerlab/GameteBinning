# /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zCodeTest_asPollinator_HQ/marker_new20191028_tuning_MQ100_with_MQ200/20191028_new_marker_set_MQ100plusManual_sorted_update.txt

marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/leaf_read_alignment_manualcuratedRef/shoremap_converted/20190819_converted_variant.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes

del_marker_finder --marker ${marker} --min-del-size 2000 --chrsizes ${sizes} -o fun

# R
# x<-read.table("del_like_20190819_converted_variant_update.bed")
# y<- x$V4-x$V3
# mean(y)
# 5511.804
# max(y)
# 351666
# min(y)
# 1000
# sum(y)
# 126584080 

# pdf("del_size_distribution.pdf")
# hist(x$V4-x$V3, breaks=1000, xlim=c(1, 50000))
# dev.off()
