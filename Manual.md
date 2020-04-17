Manual
=
This is the pipeline of how gamete binning works.

##### Prepare reads

Reads should include long reads (e.g., PacBio/Nanopore) from somatic tissue and short reads from
single cell sequencing of gamete genomes, supposing the following are ready:

* long_reads_raw.fa
* gamete_R1.fastq.gz, gamete_R2.fastq.gz


##### trim 10x barcodes


This leads to

* gamete_R1_clean.fastq.gz, gamete_R2_clean.fastq.gz


##### Genome size estimation

Count k-mers

    zcat gamete_R1_clean.fastq.gz gamete_R2_clean.fastq.gz | jellyfish count /dev/fd/0  -C -o gamete_21mer_trimmed -m 21 -t 20 -s 5G
    jellyfish histo -h 200000 -o gamete_21mer_trimmed.histo gamete_21mer_trimmed

Estimate genome size ~ 242.5 Mb
    R
    library("findGSE")
    findGSE(histo="gamete_21mer_trimmed.histo", sizek=21, outdir=".", exp_hom=200)

##### Preliminary assembly

    canu -p preasm -d canu_preasm useGrid=false genomeSize=242500000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw long_reads_raw.fa executiveThreads=20 >canu_preasm.log

or,

    flye --pacbio-raw long_reads_raw.fa --genome-size 243m --out-dir flye_preasm --threads 4

This leads to preliminary assembly

* preasm.fasta

##### Curation of assembly (with purge_haplotigs pipeline)

This leads to curated assembly

* curated_asm.fasta
