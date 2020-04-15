Gamete binning
=

Here we introduce the pipeline for gamete binning, a method based on single-cell sequencing of gamete genomes to achieve haplotype-resolved chromosome-level genome assembly of a heterozygous species.

Pre-requisite
=
Note 1: except the self-developed tools (for installation, please check INSTALL), a few publicly available tools are also needed,

* [bowtie2](https://github.com/BenLangmead/bowtie2)
* [samtools](https://github.com/samtools/)
* [bcftools](https://samtools.github.io/bcftools/)
* [SHOREmap](http://bioinfo.mpipz.mpg.de/shoremap/)
* [KMC](https://github.com/refresh-bio/KMC)
* ...

Note 2: [zlib.h](https://github.com/madler/zlib) is required by some tools (for the purpose of zipping files), please also install accordingly.
