Manual
=
This is the pipeline explaining how gamete binning works.

##### Prepare data

All data are from a single heterozygous individual of interest, including long reads (e.g., PacBio/Nanopore) from somatic tissue and short reads from
single cell sequencing of hundreds of gamete genomes. Suppose the following are ready:

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

Note, reference contig id might include "|", which is not accepted by many existing tools. Replace it with sed 's/|arrow|/_/g'.


##### Step 4. Curation of assembly (with purge_haplotigs pipeline)

Align Illumina reads to the preliminary assembly

    bowtie2 -x pre_asm_pilon.fasta -1 gamete_libx_R1_clean.fastq.gz -2 gamete_libx_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o RP_PE_sorted.bam -

Generate a coverage histogram

    purge_haplotigs readhist -b RP_PE_sorted.bam -g pre_asm_pilon.fasta -t 4

Select coverage cutoffs (./tmp_purge_haplotigs/MISC/aligned_pe.bam.histogram.csv) - here we have het-peak at 84x, hom-peak at 171x:

    lowcutoff=40   -- this is the value at the valley on the left of het-peak
    midcutoff=121  -- this is the value at the valley between het- and hom-peaks
    highcutoff=342 -- this is ~2*hom-peak

Analyse the coverage on a contig by contig basis. This script produces a contig coverage stats csv file with suspect contigs flagged for further analysis or removal.

    purge_haplotigs contigcov -i RP_PE_sorted.bam.gencov -l ${lowcutoff} -m ${midcutoff} -h ${highcutoff} -o coverage_stats.csv -j 95 -s 80

Run a BEDTools windowed coverage analysis (note, this would lead to curated.fasta and curated.haplotigs.fasta)

    ABAM=RP_PE_sorted.bam
    genome=pre_asm_pilon.fasta
    purge_haplotigs purge -g ${genome} -c coverage_stats.csv -t 16 -o curated -d -b ${ABAM} -wind_min 1000 -wind_nmax 250 -v

Re-check haplotigs (blast them with the curated genome, i.e., the one built up with selected primary contigs). If a defined haplotig is not covered by more than 50% (by primary contigs), correct it as a primary contig, and merge it with curated.fasta. (Note, visualization of the blast result with R_scripts_aux/visualize_blast_haplotig_against_purged.R -- necessary settings on paths needed). 

    db=curated.fasta
    makeblastdb -in ${db} -dbtype nucl > formatdb.log
    contigpath=/path/to/tmp_purge_haplotigs/CONTIGS/
    grep '>' curated.haplotigs.fasta | sed 's/ /\t/g' | cut -f1 | sed 's/>//' > haplotigs_as_predicted_by_purgeHaplotig.list
    while read r; do q=${r}.fasta; blastn -query ${contigpath}/${q} -db ${db} -out ${q}.oblast -outfmt 7 > blastall_blastn.log;done < ../haplotigs_as_predicted_by_purgeHaplotig.list

This leads to a version of manually curated assembly

* manually_curated.fasta

Note, although this is supposed to be a haploid assembly, it is built up with a high mixture of both haplotypes.

##### Step 5. Read alignment of pooled gamete nuclei for SNP marker definition

Index the reference sequence,

    refgenome=manually_curated.fasta
    bowtie2-build -f ${refgenome} ${refgenome} --threads 4

Align pooled gamete reads to the reference

    bowtie2 -x ${refgenome} -1 gamete_libx_R1_clean.fastq.gz -2 gamete_libx_R2_clean.fastq.gz -p 20 | samtools view -@ 20 -bS - | samtools sort -@ 20 -o gamete_ManualCurated.bam -; samtools depth gamete_ManualCurated.bam > gamete_ManualCurated.depth.txt; samtools index gamete_ManualCurated.bam

Get vcf and variants (a file of ploidy need, columns are '\t'-separated)

    echo "*	*	*	*	2" > ploidy
    bcftools mpileup -d 800 -f ${refgenome} gamete_ManualCurated.bam |bcftools call -m -v -Ov --ploidy-file ${ploidy} > gamete_ManualCurated.vcf

Convert the variant format to plain text

    SHOREmap convert --marker gamete_ManualCurated.vcf --folder shoremap_converted --indel-size 5 --min-AF 0.1 -runid 20200426
    
Filtering for allelic snps. This needs to tune quality and coverage according to the avaiable data. For allele coverage, we are looking for allelic snps, so it should be around 0.5; while the coverage on the alt allele should be around half the genome wide average. In the meanwhile, the total coverage should not be too small or too large. In my case, the avg is around 171x (and half:84x), and we can try the cutoffs below:

Note, $6 is mapping quality; $7 is coverage of alt allele, we can try with

    awk '$6>=100 && $7>=60 && $7<=140 && $7/$8>=120 && $7/$8<=280 && $8>=0.38 && $8<=0.62' /path/to/20200426_converted_variant.txt > final_snp_markers.txt

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
