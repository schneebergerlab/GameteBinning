# pb reads of currot-haplotype and orange-red-haplotye separated from rojo pasion is here - done by Jose with trio-binning method:
# /srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Triobin_scratch/analysis/seqs/q_read/Diff_kmer/ClassifyPython27/readsB 

cd /netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/zzz2019.12.13_Project_4350_CNV_Plum_Apricot_4279_resequencing_16bp_out_repeat_analysis/zTrio_binning_separation_fly_assembly/

# currot-haplotype:A_Currot

hapCU=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Triobin_scratch/analysis/seqs/q_read/Diff_kmer/ClassifyPython27/readsA/3648_A.all.subreads.bam.pb_selected.fa
bsub -m "hpc001 hpc003 hpc004 hpc005 hpc006" -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 50000 -o hpc_flye_RPhapCurrot.log -e hpc_flye_RPhapCurrot.err "flye --pacbio-raw ${hapCU} --genome-size 240m --out-dir flye_zadd_RPhapCurrot --threads 20"

# orangered-haplotype:B_OR

hapOR=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Triobin_scratch/analysis/seqs/q_read/Diff_kmer/ClassifyPython27/readsB/3648_A.all.subreads.bam.pb_selected.fa
bsub -m "hpc001 hpc003 hpc004 hpc005 hpc006" -q ioheavy -n 20 -R "span[hosts=1] rusage[mem=80000]" -M 80000 -o hpc_flye_RPhapOrangeRed.log -e hpc_flye_RPhapOrangeRed.err "flye --pacbio-raw ${hapOR} --genome-size 240m --out-dir flye_zadd_RPhapOrangeRed --threads 20"
