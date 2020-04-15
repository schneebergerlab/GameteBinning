# step 1. get PM values of contigs corrected according to contig phasing according JoinMap analysis

# Input A: where groups of contigs in 1-8.txt are collected: 

cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/

# Input B: where asPollinator works:

workpath=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_cells/z20191101_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/

# correction with the above 2 inputs:

>PM_update.log
for i in {1..8}; do s6_pollen_bed_PMflipper ${i}.txt ${workpath} ${i}; done >> PM_update.log

# where updated PM of contigs are collected: /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_cells/z20191101_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s6_PM_pollen_bed_ctgwise_JoinMap_flipped/

########################################################################################################################

# step 2. get pollen-wise genotype contigs

cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/

correctedPMpath=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_cells/z20191101_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s6_PM_pollen_bed_ctgwise_JoinMap_flipped/

# get pollen contig PM-clusters
mkdir pollen_contigs
cd pollen_contigs

for i in {1..895}; do 
   >pollen_contigs_${i}.bed
   for ctg in ${correctedPMpath}/*.txt; do
       #ll ${ctg}
       grep -P "\t"${i}"_x" ${ctg} >> pollen_contigs_${i}.bed
   done
done
#
# get contig ids of groups
resultfolder=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/
cd ${resultfolder}
mkdir grp_contig_ids
subfolder=grp_contig_ids
PMupdatedfolder=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_cells/z20191101_phasing_markers_with_correction_pollen_AB_subset895_samples_Marker_MQ100plusManualp0p38AF_scorep81_tmp_pollen_genotypes/s6_PM_pollen_bed_ctgwise_JoinMap_flipped/
for grp in {1..8}; do
    cd ${PMupdatedfolder}
    ls s6_PM_region_pollens_at_contig_*_pilon*_grp${grp}.txt | sed 's/_/\t/g' | cut -f7-8 | sed 's/\t/_/g' > ${resultfolder}/${subfolder}/grp${grp}_contigids.list
done

## for each group of contigs, find corresponding sub-pollens contigs
workfolder=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/
cd ${workfolder}
mkdir grp_wise_pollen_contigs
cd grp_wise_pollen_contigs
for grp in {1..8}; do 
    cd ${workfolder}/grp_wise_pollen_contigs
    mkdir grp${grp}_pollen_contig
    cd grp${grp}_pollen_contig
    for i in {1..895}; do    
        >pollen_contigs_${i}_grp${grp}_MMM.bed
        >pollen_contigs_${i}_grp${grp}_PPP.bed
        while read ctg; do
            grep -P ${ctg}"\t" ../../pollen_contigs/pollen_contigs_${i}.bed | awk '$6=="M" ' >> pollen_contigs_${i}_grp${grp}_MMM.bed
            grep -P ${ctg}"\t" ../../pollen_contigs/pollen_contigs_${i}.bed | awk '$6=="P" ' >> pollen_contigs_${i}_grp${grp}_PPP.bed
        done < ../../grp_contig_ids/grp${grp}_contigids.list
    done 
done

# apply below to each pollen_contigs_${i}_grp${grp}_MMM.bed, i=pollen, grp=grpid: all sub-contigs from pollen i for group grp
#bami=bamofpolleni
#samtools view -bS -L pollen_contigs_${i}_grp${grp}_MMM.bed ${bami} | samtools view -bS - | samtools sort -o pollen_contigs_${i}_grp${grp}_MMM.bam -
#samtools view -bS -L pollen_contigs_${i}_grp${grp}_PPP.bed ${bami} | samtools view -bS - | samtools sort -o pollen_contigs_${i}_grp${grp}_PPP.bam -

pollen_contigs_${i}_grp${grp}_PPP.bam

# prepare bam list of all pollens
cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/

cat -n cells_AB_subset895_list_paths.txt | awk '{printf("i=%s;ibam=%s\n", $1, $2)}' > cells_AB_subset895_list_paths_with_pollenid.txt
sed 's/sample_A\//sample_A\t/g' cells_AB_subset895_list_paths.txt | sed 's/sample_B\//sample_B\t/g' | cut -f2 | awk '{printf ("part1_%s.bam\n", $1)}' > cells_AB_subset895_list_bams.txt
paste cells_AB_subset895_list_paths_with_pollenid.txt cells_AB_subset895_list_bams.txt | sed 's/\tpart1_/\/part1_/g' > pollenid_pollenbam.list

# prepare cmds for extracting reads in bam format within each linkage group: note: as reads might not be paired anymore in bed regions, we extract 1 fastq files for each pollen MMM/PPP

workfolder=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/grp_wise_pollen_contigs/
cd ${workfolder}
bams=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/pollenid_pollenbam.list

for grp in {1..8}; do 
    cd ${workfolder}/grp${grp}_pollen_contig
    >../tmp_zgrp${grp}_MMMcluster_read_extraction_from_bams.sh  
    >../tmp_zgrp${grp}_PPPcluster_read_extraction_from_bams.sh    
    for i in {1..895}; do 
        echo "ibed="${workfolder}/grp${grp}_pollen_contig/pollen_contigs_${i}_grp${grp}_MMM.bed"; samtools view -bS -L \${ibed} \${ibam} | samtools view -bS - | samtools sort -o pollen_contigs_${i}_grp${grp}_MMM.bam -; bamToFastq -i pollen_contigs_${i}_grp${grp}_MMM.bam -fq pollen_contigs_${i}_grp${grp}_MMM_R12.fastq; gzip pollen_contigs_${i}_grp${grp}_MMM_R12.fastq" >> ../tmp_zgrp${grp}_MMMcluster_read_extraction_from_bams.sh
        echo "ibed="${workfolder}/grp${grp}_pollen_contig/pollen_contigs_${i}_grp${grp}_PPP.bed"; samtools view -bS -L \${ibed} \${ibam} | samtools view -bS - | samtools sort -o pollen_contigs_${i}_grp${grp}_PPP.bam -; bamToFastq -i pollen_contigs_${i}_grp${grp}_PPP.bam -fq pollen_contigs_${i}_grp${grp}_PPP_R12.fastq; gzip pollen_contigs_${i}_grp${grp}_PPP_R12.fastq" >> ../tmp_zgrp${grp}_PPPcluster_read_extraction_from_bams.sh
    done
    paste ${bams} ../tmp_zgrp${grp}_MMMcluster_read_extraction_from_bams.sh | sed 's/\t/; /g' > ../zgrp${grp}_MMMcluster_read_extraction_from_bams.sh
    paste ${bams} ../tmp_zgrp${grp}_PPPcluster_read_extraction_from_bams.sh | sed 's/\t/; /g' > ../zgrp${grp}_PPPcluster_read_extraction_from_bams.sh  
    rm ../tmp_zgrp${grp}_MMMcluster_read_extraction_from_bams.sh  
    rm ../tmp_zgrp${grp}_PPPcluster_read_extraction_from_bams.sh         
done
cd ${workfolder}

# run read extraction with the above cmds
workfolder=/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.10.30_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out/zPhasing_result_for_Jose_LG_analysis/phased_markers_RP_leaf_markers_0.95.nan.MARKERS_0.15_GAMETES_filtered/grp_wise_pollen_contigs/
cd ${workfolder}
for grp in {1..8}; do 
    cd ${workfolder}
    mkdir zfinal_reads_extracted_MMM_PPP_grp${grp}_allpollens
    cd zfinal_reads_extracted_MMM_PPP_grp${grp}_allpollens
    nohup ../zgrp${grp}_MMMcluster_read_extraction_from_bams.sh &>read_extraction_MMM.log&
    nohup ../zgrp${grp}_PPPcluster_read_extraction_from_bams.sh &>read_extraction_PPP.log&
done
cd ${workfolder}


































