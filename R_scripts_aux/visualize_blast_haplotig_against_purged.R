# you need to update paths below to use this script
path <- "/path/to/blast_of_haplotigs_to_curated_assembly/"
sizes   <- read.table("/path/to/contig_sizes.txt") # you need to generate a contigs size file, format: contig_name tab contig_size
primary <- read.table("/path/to/haplotigs_as_predicted_by_purgeHaplotig.list", sep="\t") # these are contigs identified as haplotigs by purge_haplotigs
#
library("scales") # transparency
#
#
pdf(paste(path, "purgedPrimary_somatic_haplotig_blast_visual_top10hit.pdf", sep=""), family="Helvetica", height=1.26772*8, width=8.26772)
par(mfrow=c(8,1), mai = c(0.6, 0.3, 0.4, 0.1)); # margin: bottom, left, top, right
#
for (ii in 1:length(primary[, 1])) # check each primary contig
{
  # current primary contig
  this_prim = as.character(primary[ii,1]) 
  # size of current primary contig
  this_size = sizes[sizes$V1==this_prim, ]
  #
  this_blas2 = read.table(paste(path, "/blast_result_primary_contigs_after_purging/", this_prim, ".fasta.oblast", sep="")) # you may need to update path setting
  this_blas  = this_blas2[this_blas2$V2!=this_prim, ]
  # unique haplotig ids
  this_hapl = unique(this_blas$V1) 
  this_matched = unique(this_blas$V2) # contigs have alignments with target haplotig 
  # check each haplotig alignment to the current primary contig
  plot(0,1, xlim=c(1, this_size[1,2]), ylim=c(-5, 10), xlab=paste(this_prim, " coord (haplotig)"), ylab="", frame.plot=F, yaxt = "n", xaxt = "n", col="white")
  lines(c(1, this_size[1,2]), c(-2.5, -2.5), lwd=30, lend = 1, col="lightgray")  # this is potential haplotig
  axis(side = 1, at = c(1, round(this_size[1,2]/2, digits = 0), this_size[1,2]))
  ## blasted segments
  interestedlen = min(length(this_matched), 10) # consider top 10 matched contigs to the above haplotig
  iicol    <- rainbow(interestedlen)
  iisizes  <- vector()
  for (jj in 1:interestedlen)
  {
    this_contig        <- as.character(this_matched[jj])
    this_hap_prim_blas <- this_blas[this_blas$V2==this_contig, ]
    iisizes            <- rbind(iisizes, paste(as.character(this_matched[jj]), sizes[sizes$V1==this_contig, 2], sep=":"))
    #
    repcol=alpha(iicol[jj], 0.6/jj)
    for (ali in c(1:length(this_hap_prim_blas$V1)))
    {
      if(this_hap_prim_blas[ali, 4]>=30) # hit >= 30 bp
      {
        alista = this_hap_prim_blas[ali, 7]
        aliend = this_hap_prim_blas[ali, 8]
        lines( c(alista, aliend),
               c(0, 0),
               col=repcol, 
               lwd=15/interestedlen*(interestedlen+1-jj),
               lend = 1)        
      }
    }  
  }
  par(xpd=TRUE)
  legend(0, 20, 
         pch    = rep(15, interestedlen),
         col    = alpha(iicol, 0.3),
         legend = iisizes,
         cex    = 0.8,
         box.col="NA",
         horiz = T,
         bg="NA"
  )
  par(xpd=FALSE)
}
##
dev.off()
