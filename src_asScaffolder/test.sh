# 

ls genetic_map_limited_markers.txt
ls zphase_contigs_linksage.txt
ls s2_genotype_contig_seq.txt
ls s2_genotype_contig_seq_del_like.txt

asCaffolder_v2 --map genetic_map_limited_markers.txt --phase zphase_contigs_linksage.txt --marker s2_genotype_contig_seq.txt --marker-del-like s2_genotype_contig_seq_del_like.txt --hap2hom 0.99 -o test0 > test0.log

# 0.99 hap2hom nearly means "do not make any change"

