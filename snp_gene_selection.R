###genetic mapping by sequencing
### downstream analysis and select the candidates genes

############################################################################
# SNPplot: plot every SNP site across on each chromosome and 
# 0/1 will be plot with red color and 1/1 with black color.

# SNPdensityPlot: plot the SNP density across each chromosome
# FilterAnn: filter the vcf file generate by snpEff package
#
############################################################################
SNPplot <- function(vcffile) { #vcffile:path of vcf file.
  library(stringr)
  vcftable <- read.table(vcffile)
  ## distinguish the genotypes 0/1 is heterozygous, 1/1 is homozygous
  ## 0/1 will be plot with red color and 1/1 with black color
  vcftable$GT <- grepl("0/1:", vcftable[,10]) 
  vcftable$col <- factor(ifelse(vcftable$GT,"yes","no"))
  chr1 <- vcftable[vcftable[,1]=="chrI",]
  chr2 <- vcftable[vcftable[,1]=="chrII",]
  chr3 <- vcftable[vcftable[,1]=="chrIII",]
  chr4 <- vcftable[vcftable[,1]=="chrIV",]
  chr5 <- vcftable[vcftable[,1]=="chrV",]
  chr6 <- vcftable[vcftable[,1]=="chrX",]
  pdf(gsub(".{3}$", "pdf", vcffile))
  par(mfrow=c(3,2),mar=c(4.1,4.1,1.1,0.1))
  plot(chr1[,2], chr1[,6], ylim=c(0,1500), col=factor(chr1$col), pch=16,
       ylab="SNP quality", xlab="", main = "ChrI")
  plot(chr2[,2], chr2[,6], ylim=c(0,1500), col=factor(chr2$col),pch=16,
       ylab="", xlab="", main = "ChrII")
  plot(chr3[,2], chr3[,6], ylim=c(0,1500), col=factor(chr3$col),pch=16,
       ylab="SNP quality", xlab="",main = "ChrIII")
  plot(chr4[,2], chr4[,6], ylim=c(0,1500), col=factor(chr4$col),pch=16,
       ylab="", xlab="",main = "ChrIV")
  plot(chr5[,2], chr5[,6], ylim=c(0,1500), col=factor(chr5$col),pch=16,
       ylab="SNP quality", xlab="", main = "ChrV")
  plot(chr6[,2], chr6[,6], ylim=c(0,1500), col=factor(chr6$col),pch=16,
       ylab="", xlab="", main = "ChrX")
  dev.off()
}


SNPdensityPlot <- function(vcffile){
  require(ggplot2)
  Cnames <- c("chr", "postion","cv","ref","mut","cv","cv2","cv3","cv4","cv5")
  F1X <- read.table(vcffile)
  names(F1X) <- Cnames
  F1X <- F1X[F1X[,6] > 200,]
  F1X$GT <- grepl("1/1:", F1X[,10]) 
  F1X <- F1X[F1X[,11],]

  SNPdensity <- ggplot(F1X) +
    geom_histogram(aes(x=postion), binwidth = 1e6) +
    facet_wrap(~ chr, ncol=2) +
    ggtitle(substr(vcffile, 59,64)) +
    xlab("Positon in genome") +
    ylab("SNP density") +
    theme_bw()
  pdf(gsub(".{10}$", "density.pdf", vcffile))
  print(SNPdensity)
  dev.off()
  
}


##############
### filter annotations file and generate a csv file contains candidates genes
###
FilterAnn <- function(vcffile, chr=c("chrI","chrII","chrIII","chrIV","chrV","chrX")){
  F1ann <- read.table(vcffile)
  F1chr <- F1ann[F1ann[,1] %in% chr,]
  F1chr$GT <- grepl("1/1", F1chr[,10])
  F1chrHomo <- F1chr[F1chr[,11],]
  F1chrHomo$NON_SYNONYMOUS_CODING <- grepl("NON_SYNONYMOUS_CODING", F1chrHomo[,8]) 
  F1chrHomo$STOP_GAINED <- grepl("STOP_GAINED", F1chrHomo[,8]) 
  F1chrHomo$FRAME_SHIFT <- grepl("FRAME_SHIFT", F1chrHomo[,8]) 
  F1chrHomo$SPLICE_SITE_ACCEPTOR <- grepl("SPLICE_SITE_ACCEPTOR", F1chrHomo[,8])
  F1chrHomo$SPLICE_SITE_DONOR <- grepl("SPLICE_SITE_DONOR", F1chrHomo[,8])
  
  F1csv <- F1chrHomo[(F1chrHomo[,12] + F1chrHomo[,13] + F1chrHomo[,14]+ F1chrHomo[,15] + F1chrHomo[,16]) >= 1,]
  for (i in 1: nrow(F1csv)) {
    test <- startsWith(unlist(strsplit(as.character(F1csv[,8][i]),",")), 
                       "NON_SYNONYMOUS_CODING")
    if(sum(test) > 0) {
      temp <- unlist(strsplit(as.character(F1csv[,8][i]),","))
      F1csv$type[i] <- temp[startsWith(temp, "NON_SYNONYMOUS_CODING")][1]
    }
    test <- startsWith(unlist(strsplit(as.character(F1csv[,8][i]),",")), 
                       "SPLICE_SITE_DONOR")
    if(sum(test) > 0) {
      temp <- unlist(strsplit(as.character(F1csv[,8][i]),","))
      F1csv$type[i] <- temp[startsWith(temp, "SPLICE_SITE_DONOR")][1]
    }
    test <- startsWith(unlist(strsplit(as.character(F1csv[,8][i]),",")), 
                       "STOP_GAINED")
    if(sum(test) > 0) {
      temp <- unlist(strsplit(as.character(F1csv[,8][i]),","))
      F1csv$type[i] <- temp[startsWith(temp, "STOP_GAINED")][1]
    }
    test <- startsWith(unlist(strsplit(as.character(F1csv[,8][i]),",")), 
                       "FRAME_SHIFT")
    if(sum(test) > 0) {
      temp <- unlist(strsplit(as.character(F1csv[,8][i]),","))
      F1csv$type[i] <- temp[startsWith(temp, "FRAME_SHIFT")][1]
    }
    test <- startsWith(unlist(strsplit(as.character(F1csv[,8][i]),",")), 
                       "SPLICE_SITE_ACCEPTOR")
    if(sum(test) > 0) {
      temp <- unlist(strsplit(as.character(F1csv[,8][i]),","))
      F1csv$type[i] <- temp[startsWith(temp, "SPLICE_SITE_ACCEPTOR")][1]
    }
    
  }
  F1csv <- F1csv[c(1:2,4:6,17)]
  names(F1csv) <- c("chr","position", "reference", "mut","quality","geneEFF" )
  write.csv(F1csv, gsub(".{3}$", "csv", vcffile), quote=F, row.names = F)
  
}


################################

## exsample
# source('snp_gene_selection.R')
# snps <- '/home/ubuntu/test/data/vcffiles/1182F1_rmdupl.raw.vcf_rmbg_homo300.vcf'
# SNPplot(snps)
# SNPdensityPlot(snps)
# 
# snpann <- '/home/ubuntu/test/data/vcffiles/1182F1_rmdupl.raw.vcf_ann_filter300.vcf'
# 
# FilterAnn(snpann)