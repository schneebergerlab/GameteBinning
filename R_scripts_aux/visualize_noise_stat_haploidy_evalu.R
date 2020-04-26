# note: you need to set up all related paths to use the script below:
#
outpath <- "/path/to/z_filter_pollen_cell_after_phasing/"
#
pdf(paste(outpath, "/check_cell_wise_PM_MP_distribution_xx.pdf", sep=""), family="Helvetica", height=12, width=8.26772)
par(mfrow=c(4,3), mai = c(0.7, 0.7, 0.5, 0.3)); # margin: bottom, left, top, right
mycols <- c("gold", "red", "blue")
#
#
########################################################################################################################
############################################# self-phased markers: YYY=578209 markers and XXX=445 Cells ################
#
# get data: the data is obtained based on markers defined with parental alleles, i.e. currot and orange red
x<-read.table("/path/to/haplotype_swaps_observed_in_cells_sorted.txt")
# for x:
# Col.1: a cell id
# Col.2: number of PM/MP observed in the cell
# Col.3: number of markers covered in the cell
########################################################################################################################
y3 <- x$V2/x$V3 # total number of PM/MP divided by total number of covered markers in a cell
hist(y3,  
     xlab="PM/MP per Covered Marker", 
     ylab="Cell Count", 
     ylim  =c(0, 400),
     breaks=seq(0,0.25, 0.01), 
     main=paste("PM/MP per Covered Marker: min~max=", min(y3), "~", round(max(y3), digits = 2), 
                "\n(Phased marker: XXX cells with YYY mkrs)", sep=""),
     cex.main=0.9,
     cex.lab =0.9,
     border = NA,
     col=mycols[3])
legend("topright",
       pch    = c(15),
       col    = "gray",
       legend = paste("Ratio<=0.05: ", sum(y3<0.05), " of ", length(y3), " cells", sep=""),
       cex    = 1,
       horiz  = F,
       box.col="NA")
write.table(file = paste(outpath, "/res_PMMP_per_Covered_Marker_LE0p05_selfphasedmkr_XXX_YYY.txt", sep=""), 
            x=x[y3<=0.05,], row.names = F, col.names = F, quote = F)
###
########################################################################################################################
# x$V2 : total number of PM/MP observed in a cell
hist(x$V2, 
     breaks = seq(0, 50000, 500), 
     xlim=c(1, 31000), 
     ylim  =c(0, 200),
     main=paste("PM/MP Count in a Cell: min~max=", min(x$V2), "~", max(x$V2), 
                "\n(Phased marker: XXX cells with YYY mkrs)", sep=""),
     cex.main=0.9,
     cex.lab =0.9,
     border = NA,
     col = mycols[3],
     xlab="PM/MP Count of a Cell", 
     ylab="Cell Count")
legend("topright",
       pch    = c(15),
       col    = "gray",
       legend = paste("Count<=2500: ", length(x[x$V2<=2500, 1]), " of ", length(x$V2), " cells", sep=""),
       cex    = 1,
       horiz  = F,
       box.col="NA")
write.table(file = paste(outpath, "/res_PMMPCount_in_Cells_LE2400_selfphasedmkr_XXX_YYY.txt", sep=""), 
            x=x[x$V2<=2500,], row.names = F, col.names = F, quote = F)
###
########################################################################################################################
plot(x$V2, 
     x$V3, 
     ylab="Number of Covered Markers of a Cell", 
     xlab="PM/MP Count of a Cell", 
     ylim  =c(0, 200000),
     xlim  =c(0, 30000),
     cex.lab=0.9,
     col=mycols[3],
     frame.plot=F)
title("PM/MP versus Covered Markers in a Cell\n XXX cells with YYY mkrs", cex.main=0.9)
###
########################################################################################################################
##
##
dev.off()

