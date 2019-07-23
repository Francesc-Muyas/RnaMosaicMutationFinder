options(warn=-1)
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-i", "--var_file", type="character", help="input read count file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-o", "--output_file", type="character", help="output_file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-cov", "--cov", type="integer", help="Minimum coverage", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()
VAR <- args$var_file
OUT <- args$output_file
COV <- args$cov

### FUNCTION
dist_SNP2 <- function(CHROM,POS,CALLS){
  distances <- abs(CALLS$POS - POS)
  distances <- distances[distances > 0]
  CLOSEST_SNP <- min(distances, na.rm = T)
  return(CLOSEST_SNP)
}

dist_SNP <- function(CALLS, SOMATIC_LIKE){
  FINAL <- ''
  NAMES <- names(table(CALLS$CHROM))
  for (chrom in NAMES){
    CALLS2 <- CALLS[CALLS$CHROM == chrom,]
    SOMATIC_LIKE2 <- SOMATIC_LIKE[SOMATIC_LIKE$CHROM == chrom,]
    
    CALLS2$DIST <- apply(CALLS2, 1, function(x)  dist_SNP2(x["CHROM"], as.numeric(x["POS"]), SOMATIC_LIKE2))
    FINAL <- rbind(FINAL, CALLS2)  
  }
  
  FINAL <- FINAL[!FINAL$CHROM == '',]
  
  return(FINAL)
}

#### LOADING FILE
SNP <- as.data.frame(fread(VAR))

#### Getting Allele Balance or Allele Frequency
SNP$AB <- SNP$ALT_COUNT/SNP$DP_HQ

#### Filtering with minimum coverage of 5
#SNP <- SNP[SNP$DP_HQ >= 5,]
SNP$FISHER <- apply(SNP,1,function(x) fisher.test(matrix(c(as.numeric(x["Ref_fwd"]), as.numeric(x["Alt_fwd"]), as.numeric(x["Ref_rev"]), as.numeric(x["Alt_rev"])), nrow = 2))[[1]])
SNP$FISHER <- p.adjust(SNP$FISHER, method = "bonferroni")

#### SOMATIC CALLING
SNP$FILTER <- (SNP$FISHER > 0.1 &
                 SNP$VariantCountBias < 5 & 
                   SNP$DP_HQ >= COV &
                   SNP$FISHER > 0.05 &
                   SNP$VariantCountBias/(SNP$ALT_COUNT + SNP$VariantCountBias) < 0.2 &
                   SNP$AB >= 0.05 & 
                   #SNP$Start_End_count/SNP$ALT_COUNT < 0.5 &
                   SNP$ALT_COUNT > 2)*1 

SNP$GT <- apply(SNP, 1, function(x) (if(as.numeric(x['FILTER']) == 1){
  "1"
  } else if (as.numeric(x['FILTER']) == 0 & as.numeric(x['DP_HQ']) >= COV & as.numeric(x['ALT_COUNT']) > 2){
    "FP"
  } else if (as.numeric(x['DP_HQ']) < COV) {
    "NA"
  } else {"."}))

#SNP[SNP$GT == "0",]$ALT <- "."

#### Printing the output
write.table(x=SNP, file=OUT, sep='\t', quote=F, row.names=F,col.names=T)

