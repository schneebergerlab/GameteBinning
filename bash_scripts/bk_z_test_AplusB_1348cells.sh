########################################################################################################################
################# step 1. preparation: check read alignments; collect sample A and B cell consensus lists ##############
workpath=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/
cd ${workpath}
sample=A
sample=B
>reads_summary.txt
while read i; do echo $i >> reads_summary.txt; grep 'reads; of these' ./sample_${sample}/$i/bowtie2.err >> reads_summary.txt; done < ./sample_${sample}/sample${sample}_cells.list
# after removing 16 "contaminated cells", we have 1345 "initial good cells", which needs to be further checking on occurrences of "P/M" swaps.
ls /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/sample_*/*/shoremap_converted/extracted_consensus_20191231_update.txt > extracted_consensus_20191231.txt
#
#
#
#
#
########################################################################################################################
################# step 2. get into working folder; prepare markers; prepare consensus list #############################
workpath=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/
cd ${workpath}
# marker
# make sure that markers on 3688 and 3826 contig is processed with contig_breaker
#
#
#
########################################################################################################################
################# step 3. select xxxx good haploid cells ################################################################
./z_test_AplusB_668cells_currot_orangered_markers.sh                            # - this is for checking selection.
zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/z_filter_pollen_cell_after_phasing/Run_check_cell_wise_potential_co_or_PMMPnoise.sh
# result with selected xxxx cells: pollen_AB_subset895_LessPMMP_samples_list.txt
#
# must confirm below is done on good cells before any asPollinator running!!!
# this is necessary for future after recreating consensus files with contig_breaker, where 3688/3826_pilon will be separated as 3688/3826_pilonsA and 3688/3826_pilonsB
# while read cc; do contig_breaker ${cc} 0; done < zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/extracted_consensus_20191231.txt
#
#
#
#
########################################################################################################################
################# step 4. phasing 1348 cells with visualization#############################################
################# application - 578209 markers further selected with MQ200, AF0.38-0.62 etc   ##############
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
date=20191231
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zCodeTest_asPollinator_HQ/marker_new20191028_tuning_MQ100_with_MQ200/20191028_new_marker_set_MQ100plusManual_sorted_update.txt
cells=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/extracted_consensus_20191231_apricot_highAlignRate.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated

bsub -q normal -R "span[hosts=1] rusage[mem=8000]" -M 8000 -o asPollinator.log -e asPollinator.err "asPollinator_v6 --marker ${marker} --pollen ${cells} -o z${date}_phasing_markers_with_correction_pollen_AB_1348_samples_Marker_MQ100plusManualp0p38AF_scorep81 --corr --ims 0.81 --size ${sizes} > z${date}_phasing_markers_with_correction_pollen_AB_1348_samples_Marker_MQ100plusManualp0p38AF_scorep81.log"

## selected 895 cells
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
date=20191231
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zCodeTest_asPollinator_HQ/marker_new20191028_tuning_MQ100_with_MQ200/20191028_new_marker_set_MQ100plusManual_sorted_update.txt
cells=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/pollen_AB_subset895_LessPMMP_samples_list_fullvarlist_20191231.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated

bsub -q normal -R "span[hosts=1] rusage[mem=8000]" -M 8000 -o asPollinator.log -e asPollinator.err "asPollinator_v6 --marker ${marker} --pollen ${cells} -o z${date}_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81 --corr --ims 0.81 --size ${sizes} > z${date}_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81.log"

# 
# Here marker set does not change
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
date=20191231
# all 1348
cd z${date}_phasing_markers_with_correction_pollen_AB_1348_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/
# subset 895
cd z${date}_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zCodeTest_asPollinator_HQ/marker_new20191028_tuning_MQ100_with_MQ200/20191028_new_marker_set_MQ100plusManual_sorted_update.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated
for i in s1_genotype*.txt; do bsub -o visualize_pollen_cluster.log -q short -R "span[hosts=1] rusage[mem=1000]" -M 2000 -m "hpc001 hpc003 hpc005 hpc006" "Rscript ../../visualize_pollen_cluster_withChrPos_v2.R ${i} ${marker} ${sizes}"; done
#
#
#
########################################################################################################################
#
#
#
#
#
#
### set up bed files for read extraction, example below:
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/z20191101_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s6_PM_pollen_bed_ctgwise/

mkdir ../s8_bed_PM
grep -P '\t196_x' *.txt | grep -P '\tP\t' | sed 's/:/\t/g' | cut -f2-9 > ../s8_bed_PM/196_x_test_P.bed
grep -P '\t196_x' *.txt | grep -P '\tM\t' | sed 's/:/\t/g' | cut -f2-9 > ../s8_bed_PM/196_x_test_M.bed

bam=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/sample_A/CCCAGTTAGCGCCTCA/part1_CCCAGTTAGCGCCTCA.bam

samtools view -bS -L 196_x_test_M.bed ${bam} | samtools view -bS - | samtools sort -o 196_x_test_M.bam -
samtools view -bS -L 196_x_test_P.bed ${bam} | samtools view -bS - | samtools sort -o 196_x_test_P.bam -





































# test unfiltered snp markers

## selected 895 cells
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
date=20191231_fullvar
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/leaf_read_alignment_manualcuratedRef/shoremap_converted/20190819_converted_variant_update.txt
cells=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/pollen_AB_subset895_LessPMMP_samples_list_fullvarlist_20191231.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated

bsub -q normal -R "span[hosts=1] rusage[mem=8000]" -M 8000 -o asPollinator.log -e asPollinator.err "asPollinator_v6 --marker ${marker} --pollen ${cells} -o z${date}_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81 --corr --ims 0.81 --size ${sizes} > z${date}_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81.log"


cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_snps/
date=20191231_fullvar
cd z${date}_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s1_genotype_pollen_seq_ctgwise/
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/leaf_read_alignment_manualcuratedRef/shoremap_converted/20190819_converted_variant_update.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes_separated
for i in s1_genotype*.txt; do bsub -o visualize_pollen_cluster.log -q short -R "span[hosts=1] rusage[mem=1000]" -M 2000 -m "hpc001 hpc003 hpc005 hpc006" "Rscript ../../visualize_pollen_cluster_withChrPos_v2.R ${i} ${marker} ${sizes}"; done






































