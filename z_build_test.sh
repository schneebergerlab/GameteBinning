cd /biodata/dep_coupland/grp_schneeberger/projects/methods/src_shq/z_asPollinator_2019/

# output phased markers using 4 of 15 good pollen samples with het-rate of < 5%

nohup ./asPollinator_v3 --marker ./marker/20190819_converted_variant_ic120_ac280_iaf0.3_aaf0.7_total250550.txt --corr --pollen ./marker/pollen4samples_list2_only4_5_11_14.txt -o z_phasing_marker --corr &>z_phasing_marker.log& 

# using phased snps to do genotyping, turn off --corr

phasedmkr=/biodata/dep_coupland/grp_schneeberger/projects/methods/src_shq/z_asPollinator_2019/z_phasing_marker_tmp_pollen_genotypes/s4_phased_markers.txt
nohup asPollinator_v2 --marker ${phasedmkr} --pollen ./marker/pollen15samples_list2.txt -o z_phased_markers_no_correction_15pollen_samples &>z_phased_markers_no_correction_15pollen_samples.log&


cd ./z_phased_markers_no_correction_15pollen_samples_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/
for i in *.txt; do bsub -o visualize_pollen_cluster.log -q short "Rscript ../../visualize_pollen_cluster_withChrPos.R ${i} ${phasedmkr}"; done 


# Graph VIZ
# dot -Tpdf -s108 s3_genotype_contig_GT_ordered_LGlike.dot > s3_genotype_contig_GT_ordered_LGlike.pdf
