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

* pre_asm.fasta

Polish the preliminary assembly using PacBio reads (using arrow/pilon), leading to the polished assembly

* pre_asm_pilon.fasta

Index it with bowtie2 for later read alignment (4 threads used)

    bowtie2-build -f pre_asm_pilon.fasta pre_asm_pilon.fasta --threads 4


##### Step 4. Curation of assembly (with purge_haplotigs pipeline)

Align Illumina reads to the reference

    bowtie2 -x pre_asm_pilon.fasta -1 gamete_libx_R1_clean.fastq.gz -2 gamete_libx_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o RP_PE_sorted.bam -

Generate a coverage histogram

    purge_haplotigs readhist -b RP_PE_sorted.bam -g pre_asm_pilon.fasta -t 4

Select coverage cutoffs (./tmp_purge_haplotigs/MISC/aligned_pe.bam.histogram.csv) - here we have het-peak at 84x, hom-peak at 171x:

    lowcutoff=40   -- this is the value at the valley on the left of het-peak
    midcutoff=121  -- this is the value at the valley between het- and hom-peaks
    highcutoff=342 -- this is ~2*hom-peak


Analyse the coverage on a contig by contig basis. This script produces a contig coverage stats csv file with suspect contigs flagged for further analysis or removal.

purge_haplotigs contigcov -i RP_PE_sorted.bam.gencov -l ${lowcutoff} -m ${midcutoff} -h ${highcutoff} -o coverage_stats.csv -j 95 -s 80

Run a BEDTools windowed coverage analysis (if generating dotplots), and 

    ABAM=RP_PE_sorted.bam
    genome=pre_asm_pilon.fasta
    purge_haplotigs purge -g ${genome} -c coverage_stats.csv -t 16 -o curated -d -b ${ABAM} -wind_min 1000 -wind_nmax 250 -v

Re-check haplotigs, 

    db=curated.fasta
    makeblastdb -in ${db} -dbtype nucl > formatdb.log
    contigpath=/path/to/tmp_purge_haplotigs/CONTIGS/
    grep '>' curated.haplotigs.fasta | sed 's/ /\t/g' | cut -f1 | sed 's/>//' > haplotigs_as_predicted_by_purgeHaplotig.list
    while read r; do q=${r}.fasta; blastn -query ${contigpath}/${q} -db ${db} -out ${q}.oblast -outfmt 7 > blastall_blastn.log;done < ../haplotigs_as_predicted_by_purgeHaplotig.list


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
