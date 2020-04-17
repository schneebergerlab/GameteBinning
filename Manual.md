Manual
=
This is the pipeline of how gamete binning works.

##### Prepare reads

Reads should include long reads (e.g., PacBio/Nanopore) from somatic tissue and short reads from
single cell sequencing of gamete genomes, supposing the following are ready:

* long_reads_raw.fa
* gamete_R1.fastq.gz, gamete_R2.fastq.gz


##### Preliminary assembly

    canu -p preasm -d canu_preasm useGrid=false genomeSize=242500000 corMhapSensitivity=high corMinCoverage=0 corOutCoverage=100 correctedErrorRate=0.105 -pacbio-raw long_reads_raw.fa executiveThreads=20 >canu_preasm.log

    flye --pacbio-raw long_reads_raw.fa --genome-size 243m --out-dir flye_preasm --threads 4

##### Curation of assembly (with purge_haplotigs pipeline)
