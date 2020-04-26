Manual
=
This is the pipeline explaining how gamete binning works.

##### Prepare data

All data are from a single heterozygous individual of interest, including long reads (e.g., PacBio/Nanopore) from somatic tissue and short reads from
single cell sequencing of hundreds of gamete genomes. Suppose the following are ready (example below):

* long_reads_raw.fa
* 4279_A_run615_SI-GA-D4_S3_L003_R1_001.fastq.gz,4279_A_run615_SI-GA-D4_S3_L003_R2_001.fastq.gz

Make softlinks (for convenience. However, 10x Genomics tools need the full name, so we keep both namings),

    ln -s 4279_A_run615_SI-GA-D4_S3_L003_R1_001.fastq.gz gamete_libx_R1.fastq.gz
    ln -s 4279_A_run615_SI-GA-D4_S3_L003_R2_001.fastq.gz gamete_libx_R2.fastq.gz

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
    
Filtering for allelic snps. This needs to tune quality and coverage according to the specific data. For allele frequency, as we are looking for allelic snps, it should be around 0.5; while the coverage on the alt allele should be around half the genome wide average. In the meanwhile, the total coverage should not be too small or too large. In my case, the avg is around 171x (and half:84x), and we can try the cutoffs below:

Note, $6 is mapping quality; $7 is coverage of alt allele, we can try with

    awk '$6>=100 && $7>=60 && $7<=140 && $7/$8>=120 && $7/$8<=280 && $8>=0.38 && $8<=0.62' /path/to/20200426_converted_variant.txt > final_snp_markers.txt

##### Step 6. 10x Genomics barcode correction

We use an artifical reference from almond at chr-level (in other cases, one can select a closely-related species with chr-level assembly): 

    wget http://getentry.ddbj.nig.ac.jp/getentry?database=na&accession_number=AP019297-AP019304&filetype=gz&limit=1000&format=fasta
    gunzip fasta_na.AP019297-AP019304.txt.gz
    sed 's/AP019297|AP019297\.1/chr1/g' fasta_na.AP019297-AP019304.txt| sed 's/AP019298|AP019298\.1/chr2/g' | sed 's/AP019299|AP019299\.1/chr3/g' | sed 's/AP019300|AP019300\.1/chr4/g' | sed 's/AP019301|AP019301\.1/chr5/g' | sed 's/AP019302|AP019302\.1/chr6/g' | sed 's/AP019303|AP019303\.1/chrX/g' | sed 's/AP019304|AP019304\.1/chrY/g' | sed 's/\./ /g' | sed 's/,//g'  | sed 's/ Prunus dulcis DNA pseudomolecule /_/g' > almond_genome.fa

Correspondingly, we prepare "a JSON file - /file_aux/contig_defs.json - describing primary contigs", and then index the genome with cellranger-dna,

    refgenome=almond_genome.fa
    cellranger-dna mkref ${refgenome} /path/to/file_aux/contig_defs.json

This will create a new reference folder of "/refdata-almond_genome/".

Correct 10x Genomics barcodes (note, if there are multiple libraries, this step needs to be done library by library, as same barcodes might be shared across libraries. However, different runs of the same library can be run together by setting option --sample=libx-run-1,libx-run-2\[,...\]),

    cellranger-dna cnv --id=4279_A_run615_cellranger --reference=/path/to/refdata-almond_genome/ --fastq=/path/to/gamete_raw_reads/ --sample=4279_A_run615_SI-GA-D4 --localcores=20 --localmem=30

Sort the bam (from the above, which is with corrected barcode information) with read name

    bam=/path/to/4279_A_run615_cellranger/outs/possorted_bam.bam
    samtools sort -n ${bam} -o RNsorted_bam.bam

Create reads ligated with corrected barcode,

    bam=/path/to/4279_A_run615_cellranger/outs/RNsorted_bam.bam
    samtools view ${bam} | T10xbam2fq - 4279_A

This leads to 

    4279_A_fqfrom10xBam_bxCorrected_R1.fq.gz
    4279_A_fqfrom10xBam_bxCorrected_R2.fq.gz

Extract individual nuclei

    R1=4279_A_fqfrom10xBam_bxCorrected_R1.fq.gz
    R2=4279_A_fqfrom10xBam_bxCorrected_R2.fq.gz
    barcode_len=16
    minimumRP=5000
    min_readpair=10000000
    asCellseparator ${barcode_len} ${R1} ${R2} ${minimumRP} ${min_readpair} cells_sep

This would report something like below,

   Info: R1 file: 4279_A_fqfrom10xBam_bxCorrected_R1.fq.gz
   Info: total number of barcodes observed in file1: 8656307
   Info: total number of reads    observed in file1: 124447727
   	among the above, good with read pairs >= 5000: 
   	number of barcodes observed in file1: 401
   	number of reads	   observed in file1: 103317664
   Info: number of effective barcodes: 401

Get list of good barcodes,

    cd /path/to/cells_sep
    awk '$2>=5000' asCellseparator_intermediate_raw_barcode_stat.txt > barcode_over_5000rpairs.list

##### step 7. Read alignment of individual gamete nuclei

Align reads of each nuclei to the curated assembly

##### Step 8. Variant calling (individual gamete nuclei)

##### Step 9. Extract allele count (individual gamete nuclei) at SNP markers

##### Step 10. Phasing SNPs within gamete genomes.

##### Step 11. Contig grouping and genetic mapping using JoinMap4.0

##### Step 12. Deletion marker definition

##### Step 13. Genetic map completing

##### Step 14. Long read separation

##### Step 15. Independent haplotype assemblies within each linkage group
