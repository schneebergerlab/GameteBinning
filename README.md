Gamete binning
=

Here we introduce the pipeline for gamete binning, a method based on single-cell sequencing of gamete genomes to achieve chromosome-level and haplotype-resolved genome assembly of a heterozygous species. This is to support the analysis in [Campoy_and_Sun_et_al_2020](https://doi.org/10.1101/2020.04.24.060046).

Pre-requisite
=
Note 1: [zlib.h](https://github.com/madler/zlib) is required by some tools (for the purpose of zipping files), please install accordingly.

Note 2: a few publicly available tools are needed,

* [bcftools](https://samtools.github.io/bcftools/)
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* [bowtie2](https://github.com/BenLangmead/bowtie2)
* [Canu v1.7](https://github.com/marbl/canu)
* [cellranger-dna](https://support.10xgenomics.com/single-cell-dna/software/pipelines/latest/what-is-cell-ranger-dna)
* [findGSE](https://github.com/schneebergerlab/findGSE)
* [Flye 2.6-g6f887ae](https://github.com/fenderglass/Flye)
* [jellyfish 2.2.10](https://github.com/gmarcais/Jellyfish)
* [KMC](https://github.com/refresh-bio/KMC)
* [purge_haplotigs](https://github.com/skingan/purge_haplotigs_multiBAM)
* [pbalign](https://github.com/PacificBiosciences/pbalign)
* [pilon](https://github.com/broadinstitute/pilon)
* [samtools](https://github.com/samtools/)
* [SHOREmap](http://bioinfo.mpipz.mpg.de/shoremap/)
* ...

Installation of self-developed tools
=

Please check [INSTALL](https://github.com/schneebergerlab/GameteBinning/blob/master/INSTALL).

Pipeline of assembly
=

Please check [Pipeline.md](https://github.com/schneebergerlab/GameteBinning/blob/master/Pipeline.md)
