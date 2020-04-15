# see snp-marker phasing at: /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis_repeat_analysis/zPhasing_cells_snps/
#
########################################################################################################################
#      
# Here it is:
#
#	1. phasing regions as deletion-like; 
#	2. separating pacbio reads; 
#	3. assembling the haplotype genome; 
#	4. scaffolding to chr-level; 
#	5. validating synteny
#
# tools: 
#
#        del_marker_finder
#        del_depth_finder
#        del_marker_genotyper
#        flye etc
#
########################################################################################################################
#
#
########################################################################################################################
# step 1. get deletion like markers according to gaps between (raw) snp variants.
#         note: marker below includes raw variants called from RP -- TODO with filtered snp markers for phasing???
#         here, also get those in the purged and manually corrected assembly but not having any variations defined.
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
marker=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/leaf_read_alignment_manualcuratedRef/shoremap_converted/20190819_converted_variant.txt
sizes=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes
#
mkdir define_del_like_regions
cd define_del_like_regions
#
del_marker_finder --marker ${marker} --min-del-size 1000 --chrsizes ${sizes} -o del_like > define_del_isize1000.log
mv del_like_isize1000_20190819_converted_variant.txt del_like_isize1000_20190819_converted_variant.bed
#
del_marker_finder --marker ${marker} --min-del-size 2000 --chrsizes ${sizes} -o del_like > define_del_isize2000.log
mv del_like_isize2000_20190819_converted_variant.txt del_like_isize2000_20190819_converted_variant.bed
#
#
R
x<-read.table("del_like_isize1000_20190819_converted_variant.bed")
y<- x$V3-x$V2
mean(y)
# 5511.804
max(y)
# 351666
min(y)
# 1000
sum(y)
# 126584080 
pdf("del_like_isize1000_20190819_converted_variant_size_distribution.pdf")
hist(y, breaks=1000, xlim=c(1, 50000), xlab="\"Del-like\" size", main="Distribution of \"del\" size (/regions without markers)")
legend("topright",
       legend = c(paste("minim: ", min(y),             " bp", sep=""), 
                  paste("maxim: ", max(y),             " bp", sep=""),        
                  paste("meann: ", round(mean(y)),     " bp", sep=""), 
                  paste("median: ",round(median(y)),   " bp", sep=""),
                  paste("<005k: ", sum(y<5000), ":",   sum(y[y<5000])/1000000,   " mb", sep=""),                                    
                  paste("<010k: ", sum(y<10000), ":",  sum(y[y<10000])/1000000,  " mb", sep=""),                  
                  paste("<020k: ", sum(y<20000), ":",  sum(y[y<20000])/1000000,  " mb", sep=""),                  
                  paste("<050k: ", sum(y<50000), ":",  sum(y[y<50000])/1000000,  " mb", sep=""), 
                  paste(">100k: ", sum(y>100000), ":", sum(y[y>100000])/1000000, " mb", sep=""),                                                    
                  paste("total: ", sum(y>0), ":",      sum(y)/1000000,           " mb", sep="")
                 ),
       horiz  = FALSE,
       border = "NA",
       bty    = "n",
       cex    = 0.8)
dev.off()
#
#
#
########################################################################################################################
# step 2. check seq depth in the above bed according to leaf short reads!
# note option -a of samtools depth must be on to used coordiate as key to match regions by del_depth_finder
# 
# 2.1 Illumina
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir depth_illumina_bam
cd depth_illumina_bam
bam=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/leaf_read_alignment_manualcuratedRef/apricot_RP_ql-trimmed_PE_sorted_ApricotManualCurated.bam
# size 1kb del-like: note -a must be on, to used coordiate as key to match regions by del_depth_finder
# bsub -o depth_bed.log -e depth_bed.err "samtools depth -a -b ../define_del_like_regions/del_like_isize1000_20190819_converted_variant.bed ${bam} > del_like_isize1000_20190819_converted_variant.bed.depth.illumina; samtoolsDepthHisto del_like_isize1000_20190819_converted_variant.bed.depth.illumina; Rscript /biodata/dep_coupland/grp_schneeberger/projects/methods/src_shq/histo_plot/histo_plot2.R del_like_isize1000_20190819_converted_variant.bed.depth.illumina.histo Apricot1kbDelLike"
# size 2kb del-like: note -a must be on, to used coordiate as key to match regions by del_depth_finder
bsub -o depth_bed.log -e depth_bed.err "samtools depth -a -b ../define_del_like_regions/del_like_isize2000_20190819_converted_variant.bed ${bam} > del_like_isize2000_20190819_converted_variant.bed.depth.illumina; samtoolsDepthHisto del_like_isize2000_20190819_converted_variant.bed.depth.illumina; Rscript /biodata/dep_coupland/grp_schneeberger/projects/methods/src_shq/histo_plot/histo_plot2.R del_like_isize2000_20190819_converted_variant.bed.depth.illumina.histo Apricot2kbDelLike"
# 2.2 PacBio (suppose bam file is ready; otherwise, align pacbio to purged version of assembly with which snps for phasing stage 1 are also defined)
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir depth_pacbio_bam
cd depth_pacbio_bam
bam=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_RP_PacbioRead_to_manually_purged_assembly/RP_pacbio_manual_purged_ref.bam
bsub -o depth_bed.log -e depth_bed.err "samtools depth -a -b ../define_del_like_regions/del_like_isize2000_20190819_converted_variant.bed ${bam} > del_like_isize2000_20190819_converted_variant.bed.depth.pacbio; samtoolsDepthHisto del_like_isize2000_20190819_converted_variant.bed.depth.pacbio; Rscript /biodata/dep_coupland/grp_schneeberger/projects/methods/src_shq/histo_plot/histo_plot2.R del_like_isize2000_20190819_converted_variant.bed.depth.pacbio.histo Apricot2kbDelLike"
# 
# 
#
######################################################################################################################## 
# step 3. check sequencing depth at del-like regions of >=2kb.
# size 2kb del-like => del_like_isize2000_20190819_converted_variant_del_like_interval_avg_depth_sorted.bed
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
cd depth_illumina_bam
region=../define_del_like_regions/del_like_isize2000_20190819_converted_variant.bed
depth=del_like_isize2000_20190819_converted_variant.bed.depth.illumina
# 3.1 illumina: hap-peak at ~220
del_depth_finder --region ${region} --depth ${depth} --hap-hom-cutoff 200 -o del_like_isize2000_20190819_converted_variant > del_depth_finder_illumina.log&
# 3.2 pacbio: hap-peak at ~33 -- caution: depth of this sequencing might be too low to trust with high confidence.
cd depth_pacbio_bam
region=../define_del_like_regions/del_like_isize2000_20190819_converted_variant.bed
depth=del_like_isize2000_20190819_converted_variant.bed.depth.pacbio
del_depth_finder --region ${region} --depth ${depth} --hap-hom-cutoff 40 -o del_like_isize2000_20190819_converted_variant > del_depth_finder_pacbio.log&
# 
# 
# 
########################################################################################################################
# step 4. get read number in the above del-like regions according to pollen samples, and normalize it as rpkm later!
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/
workpath=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/
cd ${workpath}
echo "cd ${workpath}" >4.sub_step_4_bedtools_coverage.sh

for sample in A B; do
        while read r; do printf "#${r}_${sample}\n bsub -q short -o bedtools_coverage_dellike.log -e bedtools_coverage_dellike.err \"cd ${workpath}/sample_${sample}/${r}; bedtools coverage -counts -a ../../zPhasing_cells_deletion_like/define_del_like_regions/del_like_isize2000_20190819_converted_variant.bed -b part1_${r}.bam -bed > ${r}_del_like_read_count.bed\"\n"; done < ./sample_${sample}/sample${sample}_cells.list >> 4.sub_step_4_bedtools_coverage.sh
done
#
#
#
########################################################################################################################
# step 5. genotype such del-like markers with pollen-illumina read depth and leaf-illumina read depth, and insert contigs not in linkage group into the linkage group according to phased genotypes.
#         => ./20200108_tmp_pollen_del_like_genotypes_addi/s2_genotype_contig_seq_del_like.txt
# inputs: 
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
date=20200108
leaf_depth=./depth_illumina_bam/del_like_isize2000_20190819_converted_variant_del_like_interval_avg_depth_sorted.bed
cells=../pollen_AB_subset895_LessPMMP_samples_list_fullvarlist_20191231.txt
bsub -q normal -R "span[hosts=1] rusage[mem=8000]" -M 8000 -o del_marker_genotyper.log -e del_marker_genotyper.err "del_marker_genotyper --pollen ${cells} --leaf-depth ${leaf_depth} --sample 8 --barcode 9 -o ${date} > ${date}_del_marker_genotyper.log"
#
#
#
#########################################################################################################################
# step 6. define del-like markers using asCaffolder_v2
#         => phased_s2_genotype_contig_seq_del_like.txt
#         => z_genetic_maps_updated_with_PMsimilarity_of_snp_plus_del_like_contigs/upd_map_group[1-8].txt
#            this folder includes more complete genetic maps by inserting del-like markers which will be used to anchor contigs of final assembly into chr-level.
#            CAUTION: upd_map_group[1-8].txt needs manual work: 
#                     find out contigs with left/right marker NOT next to each other; put them together at one position with minimum genetic distance to nearby markers.
#                          e.g., cat upd_map_group1.txt | cut -f1 | uniq -c | awk '$1==1' # => those showing 1 need moving. => final_manual_upd_map_group[1-8].txt
#                          s1_genotype_pollen_seq_contig_tig00003860_pilon.pdf shows a "weird" pattern.
#                     how to define a good measure?
# Input:
# zphase_contigs_linksage.txt -------- grouping   by    JoinMap
# s2_genotype_contig_seq.txt  -------- PM         by    asPollinator_v6
# s2_genotype_contig_seq_del_like.txt--del marker by    del_marker_genotyper
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
del_like=./20200108_tmp_pollen_del_like_genotypes_addi/s2_genotype_contig_seq_del_like.txt
gmap=./genetic_maps_limited_markers/genetic_map_limited_markers.txt
asCaffolder_v2 --map ${gmap} --phase zphase_contigs_linksage.txt --marker s2_genotype_contig_seq.txt --marker-del-like ${del_like} -o phased > del_phased.log
#
#
#
########################################################################################################################
# step 7. separate reads # linkage info from JoinMap # Phased snp markers from asPollinator # pb-sam file minimap2
#         ln -s /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_cells/z20191104_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s4_phased_markers.txt s4_phased_markers.txt
#         ln -s /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_RP_PacbioRead_to_manually_purged_assembly/RP_pacbio_manual_purged_ref.sam RP_pacbio_manual_purged_ref_ln.sam
#
# Input:
# s4_phased_markers.txt ------------------------ phased snps by asPollinator_v6
# phased_s2_genotype_contig_seq_del_like.txt --- phased dels by asCaffolder_v2
# zphase_contigs_linksage.txt ------------------ grouping    by JoinMap
# 
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir pbread_separation_by_pacbio_genotyper
cd pbread_separation_by_pacbio_genotyper
bsub -q normal -R "span[hosts=1] rusage[mem=3000]" -M 5000 -o pacbio_separate.log -e pacbio_separate.err  "pacbio_genotyper --sam RP_pacbio_manual_purged_ref_ln.sam --marker s4_phased_markers.txt --marker2 ../phased_s2_genotype_contig_seq_del_like.txt --phase ../zphase_contigs_linksage.txt --ims 0.9 -o rojopasionpb_separation_with_dels > rojopasionpb_intermediate_for_checkings.txt"
#
#
#
########################################################################################################################
# step 8. assembly of haplotypes.
#
# I did it with both canu and flye -- 
#   canu slow (~1 day 20 cpu) with "~100"kb higher N50 but also many small contigs; 
#   flye fast (~5 hours 20 cpu) with more mediate size contigs.
#
# step 8a1. assembly of pacbio reads with canu: 
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir assembly_20200120
cd assembly_20200120
#
path2reads=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/pbread_separation_by_pacbio_genotyper/rojopasionpb_separation_with_dels_snp_marker_separated_pbreads/
# if you want to run on hpc
#for chr in {1..8}; do
#   for sample in ${chr}.txt_PPP_pbreads.fa ${chr}.txt_MMM_pbreads.fa; do
#           cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120/
#	   bsub -m "hpc001 hpc003 hpc004 hpc005 hpc006" -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 50000 -o hpc_canu_${sample}.log -e hpc_canu_${sample}.err "/srv/netscratch/dep_coupland/grp_schneeberger/bin/canu/canu/Linux-amd64/bin//canu -p asm${sample} -d canu_${sample} useGrid=false genomeSize=30000000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw ${path2reads}/${sample} executiveThreads=20"
#   done
#done   
#
# note: cannot run all jobs on one node due to high mem requirement; start them one by one on selected dell-nodes.
# 1:dell3, 2:done,3:done,4:done,5:done,6:done,7:done,8:done
for chr in {1..8}; do
   for sample in ${chr}.txt_PPP_pbreads.fa ${chr}.txt_MMM_pbreads.fa; do
           cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120/
	   /srv/netscratch/dep_coupland/grp_schneeberger/bin/canu/canu/Linux-amd64/bin//canu -p asm${sample} -d canu_${sample} useGrid=false genomeSize=30000000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw ${path2reads}/${sample} executiveThreads=20 &>canu_chr${chr}_${sample}.log&
   done
done
#
## step 8a2. get canu N50
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120/
for i in *.txt_*_pbreads.fa; do canu_$i; /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/bin/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl asm${i}.contigs.fasta 20000000 1 > asm${i}.contigs.N50.calc.result.txt; cd ..; done
for i in */*N50*; do echo $i>>canu_N50_all.txt; cat $i>>canu_N50_all.txt; done
for i in */*PPP*N50*; do echo $i>>canu_N50_allPPP.txt; cat $i>>canu_N50_allPPP.txt; done
for i in */*MMM*N50*; do echo $i>>canu_N50_allMMM.txt; cat $i>>canu_N50_allMMM.txt; done
#
# step 8a3. canu: assemble mapped_no_lg_pb_reads.fa:57955 and unmapped_starcigar_pb_reads.fa:390453
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120/
#
path2reads=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/pbread_separation_by_pacbio_genotyper/rojopasionpb_separation_with_dels_snp_marker_separated_pbreads/
for sample in mapped_no_lg_pb_reads.fa unmapped_starcigar_pb_reads.fa; do
       cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120/
       bsub -m "hpc001 hpc003 hpc004 hpc005 hpc006" -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 50000 -o hpc_canu_${sample}.log -e hpc_canu_${sample}.err "/srv/netscratch/dep_coupland/grp_schneeberger/bin/canu/canu/Linux-amd64/bin//canu -p asm${sample} -d canu_${sample} useGrid=false genomeSize=30000000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw ${path2reads}/${sample} executiveThreads=20 stopOnReadQuality=false"       
done
#
# step 8b1. assembly of pacbio reads with Flye with assembly stat checking.
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir assembly_20200120_flye
cd assembly_20200120_flye
path2reads=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/pbread_separation_by_pacbio_genotyper/rojopasionpb_separation_with_dels_snp_marker_separated_pbreads/
for chr in 1 2 3 4 5 6 7 8; do
   for sample in ${chr}.txt_PPP_pbreads.fa ${chr}.txt_MMM_pbreads.fa; do
           cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120_flye/
	   bsub -q ioheavy -n 4 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -o hpc_flye_${sample}.log -e hpc_flye_${sample}.err "flye --pacbio-raw ${path2reads}/${sample} --genome-size 40m --out-dir flye_${sample} --threads 4"
   done
done
# step 8b2. get N50
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120_flye
for i in flye*.fa; do cd $i; /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/bin/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl assembly.fasta 20000000 1 > ${i}_N50_calc.txt; cd .. ; done
for i in */*N50*; do echo $i>>flye_N50_all.txt; cat $i>>flye_N50_all.txt; done
for i in */*PPP*N50*; do echo $i>>flye_N50_allPPP.txt; cat $i>>flye_N50_allPPP.txt; done
for i in */*MMM*N50*; do echo $i>>flye_N50_allMMM.txt; cat $i>>flye_N50_allMMM.txt; done
#
# step 8b3. flye: assemble mapped_no_lg_pb_reads.fa:57955 and unmapped_starcigar_pb_reads.fa:390453
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
cd assembly_20200120_flye
path2reads=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/pbread_separation_by_pacbio_genotyper/rojopasionpb_separation_with_dels_snp_marker_separated_pbreads/
for sample in mapped_no_lg_pb_reads.fa unmapped_starcigar_pb_reads.fa; do
       cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/assembly_20200120_flye/
       bsub -q ioheavy -n 4 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -o hpc_flye_${sample}.log -e hpc_flye_${sample}.err "flye --pacbio-raw ${path2reads}/${sample} --genome-size 40m --out-dir flye_zadd_${sample} --threads 4"
done
#
#
#
########################################################################################################################
# step 9. scaffolding contigs of purged version of assembly into chr-level.
#         total_893contigs_in_genetic_map
#
# 
## step 9a. this is scaffolded with genetic map
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir z_scaffolding_with_final_genetic_map 
cd z_scaffolding_with_final_genetic_map
fa=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.fasta
updatedgm=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_genetic_maps_updated_with_PMsimilarity_of_snp_plus_del_like_contigs
# note the i below is for labelling genetic map; different from below labelling phase group later
>scaffolder.log
for i in {1..8}; do scaffolder ${updatedgm}/final_manual_upd_map_group${i}.txt ${fa} >> scaffolder.log; done
# step 9b. check size in genetic map: 223971482 bp / 893 contigs
for i in {1..8}; do 
  >grp${i}.sizes
  cut -f1 ${updatedgm}/final_manual_upd_map_group${i}.txt | sort | uniq > gr${i}.id
  while read r; do grep $r /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/z_2019_our_Apricot_manually_curated_primary_contigs/manual_curation_our_ApricotRef/manually_curated.chrsizes >> grp${i}.sizes; done < gr${i}.id
  awk '{s+=$2} END {print s}' grp${i}.sizes
done
#
# tiny summary:
# 939 is the number of contigs in the purged assembly, totaling in size 230931811 bp. Among these 939 contigs, there are contigs showing homozygous between two parents, deletions in one of the parents, or repeats with high coverage, where we can not define a single marker; but we did it with del-like markers.
# 893 (indeed 891, as two contigs broken into 2) is the number of contigs we can find/define snp markers (at least one), totaling in size 223971482 bp.
# 653 (indeed 655, as two contigs broken into 2) is the number of contigs in the phased linkage groups (by JoinMap), totaling in size ~213 Mb.
#
#
# now scaffolding the pacbio-assembled according to the genetic-scaffolded.
#   Info: phased linkage group 1.txt matching genetically mapped group map_group2.txt; linkage files named "1.txt_MMM/PPP_pbreads.fa"
#   Info: phased linkage group 2.txt matching genetically mapped group map_group7.txt
#   Info: phased linkage group 3.txt matching genetically mapped group map_group3.txt
#   Info: phased linkage group 4.txt matching genetically mapped group map_group4.txt
#   Info: phased linkage group 5.txt matching genetically mapped group map_group6.txt
#   Info: phased linkage group 6.txt matching genetically mapped group map_group1.txt
#   Info: phased linkage group 7.txt matching genetically mapped group map_group5.txt
#   Info: phased linkage group 8.txt matching genetically mapped group map_group8.txt
#
# step 9c. prepare info for ragoo
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_scaffolding_with_final_genetic_map
#
# Caution: only when fasta files at this location can Rogoo find them!!!
for i in {1..8}; do ln -s ../assembly_20200120_flye/flye_${i}.txt_MMM_pbreads.fa/assembly.fasta ln_asm${i}.txt_MMM_pbreads.fa.contigs.fasta; done
for i in {1..8}; do ln -s ../assembly_20200120_flye/flye_${i}.txt_PPP_pbreads.fa/assembly.fasta ln_asm${i}.txt_PPP_pbreads.fa.contigs.fasta; done
#
#
# step 9d. ragoo assembly with genetic-map scaffolded reference
#
## CAUTION: genetic-map scaffolds - path setting below required by ragoo
## query=pacbio assembly for linkage group 1.txt to ref=genetic scaffolds 2
refmapfa=../scaffolded_final_manual_upd_map_group2.fa
querylab=asm1
## query=pacbio assembly for linkage group 2.txt to ref=genetic scaffolds 7
refmapfa=../scaffolded_final_manual_upd_map_group7.fa
querylab=asm2
## query=pacbio assembly for linkage group 3.txt to ref=genetic scaffolds 3
refmapfa=../scaffolded_final_manual_upd_map_group3.fa
querylab=asm3
## query=pacbio assembly for linkage group 4.txt to ref=genetic scaffolds 4
refmapfa=../scaffolded_final_manual_upd_map_group4.fa
querylab=asm4
## query=pacbio assembly for linkage group 5.txt to ref=genetic scaffolds 6
refmapfa=../scaffolded_final_manual_upd_map_group6.fa
querylab=asm5
## query=pacbio assembly for linkage group 6.txt to ref=genetic scaffolds 1
refmapfa=../scaffolded_final_manual_upd_map_group1.fa
querylab=asm6
## query=pacbio assembly for linkage group 7.txt to ref=genetic scaffolds 5
refmapfa=../scaffolded_final_manual_upd_map_group5.fa
querylab=asm7
## query=pacbio assembly for linkage group 8.txt to ref=genetic scaffolds 8
refmapfa=../scaffolded_final_manual_upd_map_group8.fa
querylab=asm8
#
for PM in PPP MMM; do  
    mkdir ${querylab}_scaffolding_${PM}
    cd ${querylab}_scaffolding_${PM}
    queryfa=ln_${querylab}.txt_${PM}_pbreads.fa.contigs.fasta # path required by ragoo
    bsub -o ragoo.log -e ragoo.err -q short -R "span[hosts=1] rusage[mem=15000]" -M 15000 "ragoo.py ../${queryfa} ${refmapfa} -C"
    cd ..
done
#
# step 9e.scaffold N50 etc: note N50 is not meaningful here as we align contigs to "reference" (from genetic map scaffolds)
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_scaffolding_with_final_genetic_map
#
path=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_scaffolding_with_final_genetic_map
for i in {1..8}; do 
    cd ${path}/asm${i}_scaffolding_MMM/ragoo_output/
    perl /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/bin/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl ragoo.fasta 30000000 1 > asm${i}_MMM_stat.txt
    cp asm${i}_MMM_stat.txt ../../
    cd ${path}/asm${i}_scaffolding_PPP/ragoo_output/
    perl /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/bin/GitHub-schneeberger-group-jac/toolbox/Analysis/Assembly/calc_CN50.pl ragoo.fasta 30000000 1 > asm${i}_PPP_stat.txt
    cp asm${i}_PPP_stat.txt ../../    
done
# artifical combination of all P and M
cat asm*_PPP_stat.txt | grep 'ContigLength' | sed 's/:\t\t/\t/g' > am1-8_contigLength_allPPP.txt 
awk '{s+=$2} END {print s}' am1-8_contigLength_allPPP.txt # => 215449629 bp
cat asm*_MMM_stat.txt | grep 'ContigLength' | sed 's/:\t\t/\t/g' > am1-8_contigLength_allMMM.txt 
awk '{s+=$2} END {print s}' am1-8_contigLength_allMMM.txt # => 214358789 bp
# real combination of all currot and orange red (combined according to kmer checking; in reality, we may not have this info)
cat asm1_PPP_stat.txt asm2_MMM_stat.txt asm3_PPP_stat.txt asm4_PPP_stat.txt asm5_MMM_stat.txt asm6_PPP_stat.txt asm7_MMM_stat.txt asm8_MMM_stat.txt | grep 'ContigLength' | sed 's/:\t\t/\t/g' > am1-8_contigLength_allOrangered_real.txt
awk '{s+=$2} END {print s}' am1-8_contigLength_allOrangered_real.txt # => 214550436
cat asm1_MMM_stat.txt asm2_PPP_stat.txt asm3_MMM_stat.txt asm4_MMM_stat.txt asm5_PPP_stat.txt asm6_MMM_stat.txt asm7_PPP_stat.txt asm8_PPP_stat.txt | grep 'ContigLength' | sed 's/:\t\t/\t/g' > am1-8_contigLength_allCurrot_real.txt
awk '{s+=$2} END {print s}' am1-8_contigLength_allCurrot_real.txt    # => 215257982
#
cat asm*_*_stat.txt | grep 'N50' | grep -v 'CN50' | sed 's/:\t\t/\t/g' | awk '{s+=$2} END {print s/16}' # => 2.66113e+07
# total N: 186200 in 16 fasta files.
>z_total_in_ragoo_fasta
for i in asm*_scaffolding_*; do 
   echo ${i} >>z_total_in_ragoo_fasta
   fasta_letter_counter ${i}/ragoo_output/ragoo.fasta N | grep 'Total number of N is: ' >> z_total_in_ragoo_fasta
done
#
#
# copy and rename ragoo fasta
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_scaffolding_with_final_genetic_map
path=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_scaffolding_with_final_genetic_map
for i in {1..8}; do 
    cd ${path}/asm${i}_scaffolding_MMM/ragoo_output/
    cp ragoo.fasta ragoo_MMM_asm${i}.fasta
    ll *.fasta
    cd ${path}/asm${i}_scaffolding_PPP/ragoo_output/
    cp ragoo.fasta ragoo_PPP_asm${i}.fasta
    ll *.fasta
done
#
# assembly size comparison => final_assemblysize_geneticsize_comparison.R # why flye lost 8 Mb compared with genetic map scaffolded? in genetic map, del-like regions are haplotypes, but doubled (they are haplotigs in the purged version of initial assembly, increasing genome size).
#
#
#
########################################################################################################################
# step 10. synteny validation
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/
mkdir z_synteny_validation_with_ragoo_fasta 
cd z_synteny_validation_with_ragoo_fasta
#
## step 10a. get new assembly for P and M
mkdir paternal_assembly
cp ../z_scaffolding_with_final_genetic_map/asm*_scaffolding_PPP/ragoo_output/ragoo_PPP_asm*.fasta ./paternal_assembly
cd ./paternal_assembly
sed -i 's/_RaGOO/_RaGOO_asm1/g' ragoo_PPP_asm1.fasta
sed -i 's/_RaGOO/_RaGOO_asm2/g' ragoo_PPP_asm2.fasta
sed -i 's/_RaGOO/_RaGOO_asm3/g' ragoo_PPP_asm3.fasta
sed -i 's/_RaGOO/_RaGOO_asm4/g' ragoo_PPP_asm4.fasta
sed -i 's/_RaGOO/_RaGOO_asm5/g' ragoo_PPP_asm5.fasta
sed -i 's/_RaGOO/_RaGOO_asm6/g' ragoo_PPP_asm6.fasta
sed -i 's/_RaGOO/_RaGOO_asm7/g' ragoo_PPP_asm7.fasta
sed -i 's/_RaGOO/_RaGOO_asm8/g' ragoo_PPP_asm8.fasta
cat *.fasta > ragoo_PPP_asm1to8.fa
fasta_length ragoo_PPP_asm1to8.fa > ragoo_PPP_asm1to8_scaffoldsizes.txt
cd ..
#
mkdir maternal_assembly
cp /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zAssembly_test_kmer_separated_with_pb_alignment/assembly_scaffolding/asm*_scaffolding_MMM/ragoo_output/ragoo_MMM_asm*.fasta ./maternal_assembly
cd ./maternal_assembly
sed -i 's/_RaGOO/_RaGOO_asm1/g' ragoo_MMM_asm1.fasta
sed -i 's/_RaGOO/_RaGOO_asm2/g' ragoo_MMM_asm2.fasta
sed -i 's/_RaGOO/_RaGOO_asm3/g' ragoo_MMM_asm3.fasta
sed -i 's/_RaGOO/_RaGOO_asm4/g' ragoo_MMM_asm4.fasta
sed -i 's/_RaGOO/_RaGOO_asm5/g' ragoo_MMM_asm5.fasta
sed -i 's/_RaGOO/_RaGOO_asm6/g' ragoo_MMM_asm6.fasta
sed -i 's/_RaGOO/_RaGOO_asm7/g' ragoo_MMM_asm7.fasta
sed -i 's/_RaGOO/_RaGOO_asm8/g' ragoo_MMM_asm8.fasta
cat *.fasta > ragoo_MMM_asm1to8.fa
fasta_length ragoo_MMM_asm1to8.fa > ragoo_MMM_asm1to8_scaffoldsizes.txt
#
#### check which contigs belong to which chromosomes/linkage groups according to closely related species ############
#
# step 10b. index reference assemblies of closely related species
#
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/2014_P_mume_Japanese_Apricot/GCF_000346735.1_P.mume_V1.0/
bsub -o minimap2_index.log -e minimap2_index.err "minimap2 -d GCF_000346735.1_P.mume_V1.0_genomic.mmi GCF_000346735.1_P.mume_V1.0_genomic.fa"
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/2017_P_persica_Peach/GCF_000346465.2_Prunus_persica_NCBIv2/
bsub -o minimap2_index.log -e minimap2_index.err "minimap2 -d GCF_000346465.2_Prunus_persica_NCBIv2_genomic.mmi GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fa"
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/2019_P_dulcis_Almond/almondReference/
bsub -o minimap2_index.log -e minimap2_index.err "minimap2 -d fasta_na.AP019297-AP019304.mmi fasta_na.AP019297-AP019304.fa"
#
# step 10c. align manually curated assembly to the closely related species
#
## REF 1: 2014_P_mume_Japanese_Apricot - the most closely related to our apricot
##
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta
mkdir PacBio_newRP_against_2014_P_mume_Japanese_Apricot
cd PacBio_newRP_against_2014_P_mume_Japanese_Apricot
#
japricotgenome=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/2014_P_mume_Japanese_Apricot/GCF_000346735.1_P.mume_V1.0/GCF_000346735.1_P.mume_V1.0_genomic.fa
apricotgenome=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta/paternal_assembly/ragoo_PPP_asm1to8.fa
ln -s ${japricotgenome} ln_japricotgenome
ln -s ${apricotgenome} ln_apricotgenome
bsub -o ragoo.log -e ragoo.err -q multicore20 -n 2 -R "span[hosts=1] rusage[mem=10000]" -M 10000 "ragoo.py -t 2 ln_apricotgenome ln_japricotgenome"
##
## REF 2: 2017_P_persica_Peach
##
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta
mkdir PacBio_newRP_against_2017_P_persica_Peach
cd PacBio_newRP_against_2017_P_persica_Peach
#
peachgenome=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/2017_P_persica_Peach/GCF_000346465.2_Prunus_persica_NCBIv2/GCF_000346465.2_Prunus_persica_NCBIv2_genomic.fa
apricotgenome=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta/paternal_assembly/ragoo_PPP_asm1to8.fa
ln -s ${peachgenome} ln_peachgenome
ln -s ${apricotgenome} ln_apricotgenome
bsub -o ragoo.log -e ragoo.err -q multicore20 -n 2 -R "span[hosts=1] rusage[mem=10000]" -M 10000 "ragoo.py -t 2 ln_apricotgenome ln_peachgenome"
##
## REF 3: 2019_P_dulcis_Almond
##
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta
mkdir PacBio_newRP_against_2019_P_dulcis_Almond
cd PacBio_newRP_against_2019_P_dulcis_Almond
#
almondgenome=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Pollen-b/z_related_species_reference_based/2019_P_dulcis_Almond/almondReference/fasta_na.AP019297-AP019304.fa
apricotgenome=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta/paternal_assembly/ragoo_PPP_asm1to8.fa
ln -s ${almondgenome} ln_almondgenome
ln -s ${apricotgenome} ln_apricotgenome
bsub -o ragoo.log -e ragoo.err -q multicore20 -n 2 -R "span[hosts=1] rusage[mem=10000]" -M 10000 "ragoo.py -t 2 ln_apricotgenome ln_almondgenome"
#
#
# step 10d. visualize synteny: check visualize_alignment.R: final_manual_upd_map_group2_RaGOO_asm1  25105107    5048622 5206625 +   NC_024128.1_LG3 24358521    4783608 4943221
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zPhasing_cells_deletion_like/z_synteny_validation_with_ragoo_fasta/PacBio_newRP_against_2014_P_mume_Japanese_Apricot/ragoo_output/

cut -f1,3,4,6,8,9,5 contigs_against_ref.paf | grep 'final_manual_upd_map_group' > data_test










