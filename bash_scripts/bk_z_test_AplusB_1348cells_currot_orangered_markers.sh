########################################################################################################################
######################################## step 1. Collect sample A and B cell lists #####################################
#
#
# TODO: revise and rerun after 20191231
#
#
########################################################################################################################
####################### step 2. phase and cluster '1000' A+B cells from 10x single cell sequencing #####################
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
# marker
#
#
# without correction= --corr off: 
# new marker sets defined with parental Currot, Orange Red and Rojo Pasion; 3688 contig separated as two contigs with contig_breaker for both consensus file and marker file
# 1. 1348 initial good cells
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzHequan_assembly_pipeline/Step_3_build_benchmark_2/marker_creation/new_marker/rojoPasion_20190819_converted_variant_ic120_ac280_iaf0.3_aaf0.7_total250550_separated.txt
date=20191231
cells=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/extracted_consensus_20191231_apricot_highAlignRate.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated

bsub -q normal -R "span[hosts=1] rusage[mem=8000]" -M 8000 -o asPollinator.log -e asPollinator.err "asPollinator_v6 --marker ${marker} --pollen ${cells} -o z${date}_parental_currot_orangered_markers_no_correction_pollen_AB_1348_samples_scorep81 --ims 0.81 --size ${sizes} > z${date}_parental_currot_orangered_markers_no_correction_pollen_AB_1348_samples_scorep81.log"

# selected 895 "good" cells
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzHequan_assembly_pipeline/Step_3_build_benchmark_2/marker_creation/new_marker/rojoPasion_20190819_converted_variant_ic120_ac280_iaf0.3_aaf0.7_total250550_separated.txt
date=20191231
cells=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/pollen_AB_subset895_LessPMMP_samples_list_fullvarlist_20191231.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated
bsub -q normal -R "span[hosts=1] rusage[mem=8000]" -M 8000 -o asPollinator.log -e asPollinator.err "asPollinator_v6 --marker ${marker} --pollen ${cells} -o z${date}_parental_currot_orangered_markers_no_correction_pollen_AB_subset895_samples_scorep81 --ims 0.81 --size ${sizes} > z${date}_parental_currot_orangered_markers_no_correction_pollen_AB_subset895_samples_scorep81.log"

# here marker set does not change!
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
date=20191231
# all 1348
cd z${date}_parental_currot_orangered_markers_no_correction_pollen_AB_1348_samples_scorep81_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/
# subset 895
cd z${date}_parental_currot_orangered_markers_no_correction_pollen_AB_subset895_samples_scorep81_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzHequan_assembly_pipeline/Step_3_build_benchmark_2/marker_creation/new_marker/rojoPasion_20190819_converted_variant_ic120_ac280_iaf0.3_aaf0.7_total250550_separated.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated
#
for i in s1_genotype*.txt; do bsub -o visualize_pollen_cluster.log -q short -R "span[hosts=1] rusage[mem=1000]" -M 2000 -m "hpc001 hpc003 hpc005 hpc006" "Rscript ../../visualize_pollen_cluster_withChrPos_v2.R ${i} ${marker} ${sizes}"; done
