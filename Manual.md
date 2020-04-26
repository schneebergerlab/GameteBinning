Manual
=
This is the pipeline explaining how gamete binning works.

##### Prepare reads

Reads should include long reads (e.g., PacBio/Nanopore) from somatic tissue and short reads from
single cell sequencing of hundreds of gamete genomes, supposing the following are ready:

* long_reads_raw.fa
* gamete_libx_R1.fastq.gz, gamete_libx_R2.fastq.gz


##### Step 1. Trim reads

Trim 16 bp barcodes off R1's (10x Genomics library setting)

    T10X_barcode_trimmer gamete_libx_R1.fastq.gz gamete_libx_R2.fastq.gz

This leads to

* gamete_libx_R1_clean.fastq.gz, gamete_libx_R2_clean.fastq.gz


##### Step 2. Genome size estimation with cleaned reads

Count k-mers

    zcat gamete_libx_R1_clean.fastq.gz gamete_libx_R2_clean.fastq.gz | jellyfish count /dev/fd/0  -C -o gamete_21mer_trimmed -m 21 -t 20 -s 5G
    jellyfish histo -h 200000 -o gamete_21mer_trimmed.histo gamete_21mer_trimmed

Estimate genome size ~ 242.5 Mb

    R
    library("findGSE")
    findGSE(histo="gamete_21mer_trimmed.histo", sizek=21, outdir=".", exp_hom=200)

##### Step 3. Preliminary assembly

    canu -p preasm -d canu_preasm useGrid=false genomeSize=242500000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw long_reads_raw.fa executiveThreads=20 >canu_preasm.log

or,

    flye --pacbio-raw long_reads_raw.fa --genome-size 243m --out-dir flye_preasm --threads 4

This leads to preliminary assembly

* preasm.fasta

##### Step 4. Curation of assembly (with purge_haplotigs pipeline)



This leads to curated assembly

* curated_asm.fasta


##### Step 5. Read alignment of pooled gamete nuclei for SNP marker definition

##### Step 6. Read alignment of individual gamete nuclei

Extract individual nuclei

Align reads of each nuclei to the curated assembly

##### Step 7. Variant calling (individual gamete nuclei)

##### Step 8. Extract allele count (individual gamete nuclei) at SNP markers

##### Step 9. Phasing SNPs within gamete genomes.

##### Step 10. Contig grouping and genetic mapping using JoinMap4.0

##### Step 11. Deletion marker definition

##### Step 12. Genetic map completing

##### Step 13. Long read separation

##### Step 14. Independent haplotype assemblies within each linkage group
