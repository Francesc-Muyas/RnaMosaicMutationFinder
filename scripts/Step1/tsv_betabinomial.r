options(warn=-1)
options(scipen=999)

suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(data.table))
#suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(bbmle))
suppressPackageStartupMessages(library(survcomp)) # Bioconductor 

# Variant calling functions
# FUNCTION
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

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

MODEL_BETABIN_DP_HQ <- function(DATA){
  result <- tryCatch({mle2(ALT_COUNT~dbetabinom.ab(size=DP_HQ,shape1,shape2),
                           data=DATA,
                           method="Nelder-Mead",
                           skip.hessian=TRUE,
                           start=list(shape1=1,shape2=round(mean(DATA$DP_HQ))),
                           control=list(maxit=1000))}, 
                     error = function(e) {(estBetaParams(mean(DATA$ALT_COUNT/DATA$DP_HQ, na.rm = T),var(DATA$ALT_COUNT/DATA$DP_HQ, na.rm = T)))})
  PARAM1 <- ifelse(is.null(coef(result)[[1]]),result[[1]], coef(result)[[1]])
  PARAM2 <- ifelse(is.null(coef(result)[[2]]),result[[2]], coef(result)[[2]])
  
  return (c(PARAM1, PARAM2))  
} 

P_VAL <- function(SNP, a, b){
  p_val <- pzoibetabinom.ab(SNP$ALT_COUNT-1, SNP$DP_HQ, a, b, lower.tail=F)
  p_val[p_val < 0] <- 0
  return(p_val)
}

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-t", "--tumor_file", type="character", help="Tumor - Read count input file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-n", "--normal_file", type="character", help="Normal - Read count input file", nargs=1, default = NULL)
parser$add_argument("-tr", "--train_file", type="character", help="File to trian error rate distribution [optional]", nargs=1, default = NULL)
parser$add_argument("-o", "--out_file", type="character", help="output_file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-b", "--bed", type="character", help="Bed file with positions to ignore for modelling", nargs=1, default = NULL)

# reading parameters
args <- parser$parse_args()

TD <- args$tumor_file

ND <- args$normal_file

TR <- args$train_file
  
bed <- args$bed 

OUTF <- args$out_file

print (paste("Normal file:", ND, sep=" "))
print (paste("Tumor file:", TD, sep=" "))
print (paste("Out file:", OUTF))


# Training set
if (is.null(TR)){
  CONTROL <- as.data.frame(fread(TD))
  CONTROL$AB <- CONTROL$ALT_COUNT/CONTROL$DP_HQ
  CONTROL <- CONTROL[CONTROL$AB < 0.1,] # Focusing only to non-germline sites
                     
  if (nrow(CONTROL) > 500000) {
    CONTROL <- CONTROL[sample(nrow(CONTROL), 500000),]
  }

    if (!is.null(bed)){
  
    BED <- as.data.frame(fread(bed,header = F))
    colnames(BED)[1] <- 'CHROM'
    colnames(BED)[2] <- 'POS'
    
    CONTROL <- CONTROL[CONTROL$AB < 0.1 & CONTROL$DP_HQ > 10 & !paste(CONTROL$CHROM,CONTROL$POS,sep = ',') %in% paste(BED$CHROM,BED$POS,sep = ','),]
    
  } else {
    CONTROL <- CONTROL[CONTROL$AB < 0.1 & CONTROL$DP_HQ > 10 & CONTROL$ALT_COUNT <= quantile(CONTROL$ALT_COUNT,0.99),]
  }
} else {
  CONTROL <- as.data.frame(fread(TR))
  CONTROL$AB <- CONTROL$ALT_COUNT/CONTROL$DP_HQ
  CONTROL <- CONTROL[CONTROL$AB < 0.1,]
}

#CONTROL <- CALLS[CALLS$AB < 0.1 & CALLS$DP_HQ > 0 & CALLS$AB > 0,] 
#CONTROL <- CONTROL[CONTROL$AB < 0.1 & CONTROL$DP_HQ > 10 & CONTROL$AB > 0,] 
#CONTROL <- CONTROL[CONTROL$AB < 0.1 & CONTROL$DP_HQ > 10,] 

## NUCLEOTIDE CHANGES
TG <- c("T>G","A>C")
TA <- c("T>A","A>T")
TC <- c("T>C","A>G")
GT <- c("G>T","C>A")
GA <- c("G>A","C>T")
GC <- c("G>C","C>G")


###### SPLIT THEM IN BASE QUALITY RANGES
# Groups of qualities

# # ONLY INCLUDING VARIANTS 
# CONTROL_TG <- CONTROL[paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  TG & CONTROL$VARIANT_TYPE == "SNP",]
# CONTROL_TA <- CONTROL[paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  TA & CONTROL$VARIANT_TYPE == "SNP",]
# CONTROL_TC <- CONTROL[paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  TC & CONTROL$VARIANT_TYPE == "SNP",]
# CONTROL_GT <- CONTROL[paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  GT & CONTROL$VARIANT_TYPE == "SNP",]
# CONTROL_GA <- CONTROL[paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  GA & CONTROL$VARIANT_TYPE == "SNP",]
# CONTROL_GC <- CONTROL[paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  GC & CONTROL$VARIANT_TYPE == "SNP",]

# # ONLY INCLUDING VARIANTS BUT INDELS SEPARATED
CONTROL_TG <- CONTROL[(paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  TG & CONTROL$VARIANT_TYPE == "SNP")
                      | (CONTROL$REF %in% c('T','A') & CONTROL$ALT_COUNT == 0),]

CONTROL_TA <- CONTROL[(paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  TA & CONTROL$VARIANT_TYPE == "SNP")
                      | (CONTROL$REF %in% c('T','A') & CONTROL$ALT_COUNT == 0),]

CONTROL_TC <- CONTROL[(paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  TC & CONTROL$VARIANT_TYPE == "SNP")
                      | (CONTROL$REF %in% c('T','A') & CONTROL$ALT_COUNT == 0),]

CONTROL_GT <- CONTROL[(paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  GT & CONTROL$VARIANT_TYPE == "SNP")
                      | (CONTROL$REF %in% c('G','C') & CONTROL$ALT_COUNT == 0),]

CONTROL_GA <- CONTROL[(paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  GA & CONTROL$VARIANT_TYPE == "SNP")
                      | (CONTROL$REF %in% c('G','C') & CONTROL$ALT_COUNT == 0),]

CONTROL_GC <- CONTROL[(paste(CONTROL$REF,CONTROL$ALT,sep='>') %in%  GC & CONTROL$VARIANT_TYPE == "SNP")
                      | (CONTROL$REF %in% c('G','C') & CONTROL$ALT_COUNT == 0),]

## Estimating parameters for each indels
CONTROL_Indel_A <- CONTROL[(CONTROL$REF %in% c('T','A'))
                           | (CONTROL$REF %in% c('T','A') & nchar(CONTROL$ALT) > 1)
                           | (CONTROL$ALT %in% c('T','A') & nchar(CONTROL$REF) > 1),]

CONTROL_Indel_G <- CONTROL[(CONTROL$REF %in% c('G','C'))
                           | (CONTROL$REF %in% c('G','C') & nchar(CONTROL$ALT) > 1)
                           | (CONTROL$ALT %in% c('G','C') & nchar(CONTROL$REF) > 1),]



## Estimating parameters for each nucleotide change when we only have 1 duplicate per barcode group
FIT_DP1.TG <- MODEL_BETABIN_DP_HQ(CONTROL_TG)

FIT_DP1.TA <- MODEL_BETABIN_DP_HQ(CONTROL_TA)

FIT_DP1.TC <- MODEL_BETABIN_DP_HQ(CONTROL_TC)

FIT_DP1.GT <- MODEL_BETABIN_DP_HQ(CONTROL_GT)

FIT_DP1.GA <- MODEL_BETABIN_DP_HQ(CONTROL_GA)

FIT_DP1.GC <- MODEL_BETABIN_DP_HQ(CONTROL_GC)

FIT_DP1.Indel_A <- MODEL_BETABIN_DP_HQ(CONTROL_Indel_A)

FIT_DP1.Indel_G <- MODEL_BETABIN_DP_HQ(CONTROL_Indel_G)


# Getting parameters
a_DP1.GA <- (FIT_DP1.GA)[[1]]
b_DP1.GA <- (FIT_DP1.GA)[[2]]

a_DP1.GC <- (FIT_DP1.GC)[[1]]
b_DP1.GC <- (FIT_DP1.GC)[[2]]

a_DP1.GT <- (FIT_DP1.GT)[[1]]
b_DP1.GT <- (FIT_DP1.GT)[[2]]

a_DP1.TA <- (FIT_DP1.TA)[[1]]
b_DP1.TA <- (FIT_DP1.TA)[[2]]

a_DP1.TC <- (FIT_DP1.TC)[[1]]
b_DP1.TC <- (FIT_DP1.TC)[[2]]

a_DP1.TG <- (FIT_DP1.TG)[[1]]
b_DP1.TG <- (FIT_DP1.TG)[[2]]

a_DP1.Indel_A <- (FIT_DP1.Indel_A)[[1]]
b_DP1.Indel_A <- (FIT_DP1.Indel_A)[[2]]

a_DP1.Indel_G <- (FIT_DP1.Indel_G)[[1]]
b_DP1.Indel_G <- (FIT_DP1.Indel_G)[[2]]

rm(CONTROL)
rm(CONTROL_GA)
rm(CONTROL_GC)
rm(CONTROL_GT)
rm(CONTROL_TA)
rm(CONTROL_TC)
rm(CONTROL_TG)
rm(CONTROL_Indel_A)
rm(CONTROL_Indel_G)

############
#### SNP Calling
############
CALLS <- as.data.frame(fread(TD))
CALLS$AB <- CALLS$ALT_COUNT/CALLS$DP_HQ

# Groups of nt changes
CALLS_TG <- CALLS[paste(CALLS$REF,CALLS$ALT,sep='>') %in%  TG,]
CALLS_TA <- CALLS[paste(CALLS$REF,CALLS$ALT,sep='>') %in%  TA,]
CALLS_TC <- CALLS[paste(CALLS$REF,CALLS$ALT,sep='>') %in%  TC,]
CALLS_GT <- CALLS[paste(CALLS$REF,CALLS$ALT,sep='>') %in%  GT,]
CALLS_GA <- CALLS[paste(CALLS$REF,CALLS$ALT,sep='>') %in%  GA,]
CALLS_GC <- CALLS[paste(CALLS$REF,CALLS$ALT,sep='>') %in%  GC,]

#### FOR GA group
CALLS_GA$P_VAL <- P_VAL(CALLS_GA, a_DP1.GA, b_DP1.GA)
CALLS_GA$P_VAL_adj <- p.adjust(CALLS_GA$P_VAL, method = "bonferroni")

#### FOR GC group
CALLS_GC$P_VAL <- P_VAL(CALLS_GC, a_DP1.GC, b_DP1.GC)
CALLS_GC$P_VAL_adj <- p.adjust(CALLS_GC$P_VAL, method = "bonferroni")

#### FOR GT group
CALLS_GT$P_VAL <- P_VAL(CALLS_GT, a_DP1.GT, b_DP1.GT)
CALLS_GT$P_VAL_adj <- p.adjust(CALLS_GT$P_VAL, method = "bonferroni")

#### FOR TA group
CALLS_TA$P_VAL <- P_VAL(CALLS_TA, a_DP1.TA, b_DP1.TA)
CALLS_TA$P_VAL_adj <- p.adjust(CALLS_TA$P_VAL, method = "bonferroni")

#### FOR TC group
CALLS_TC$P_VAL <- P_VAL(CALLS_TC, a_DP1.TC, b_DP1.TC)
CALLS_TC$P_VAL_adj <- p.adjust(CALLS_TC$P_VAL, method = "bonferroni")

#### FOR TG group
CALLS_TG$P_VAL <- P_VAL(CALLS_TG, a_DP1.TG, b_DP1.TG)
CALLS_TG$P_VAL_adj <- p.adjust(CALLS_TG$P_VAL, method = "bonferroni")

# Merge SNP calls
SNP <- rbind(CALLS_GA, CALLS_GC, CALLS_GT, CALLS_TA, CALLS_TC, CALLS_TG)


############
#### Indel Calling
############

Indel_A <- CALLS[(CALLS$REF %in% c('T','A') & nchar(CALLS$ALT) > 1 & CALLS$ALT_COUNT > 0) |
                   (CALLS$ALT %in% c('T','A') & nchar(CALLS$REF) > 1 & CALLS$ALT_COUNT > 0),]

Indel_G <- CALLS[(CALLS$REF %in% c('G','C') & nchar(CALLS$ALT) > 1 & CALLS$ALT_COUNT > 0) |
                   (CALLS$ALT %in% c('G','C') & nchar(CALLS$REF) > 1 & CALLS$ALT_COUNT > 0),]

# Indel_A
if ( nrow(Indel_A) > 0 ){
  Indel_A$P_VAL <- P_VAL(Indel_A, a_DP1.Indel_A, b_DP1.Indel_A)
  Indel_A$P_VAL_adj <- p.adjust(Indel_A$P_VAL, method = "bonferroni")
}

# Indel_G
if ( nrow(Indel_G) > 0 ){
  Indel_G$P_VAL <- P_VAL(Indel_G, a_DP1.Indel_G, b_DP1.Indel_G)
  Indel_G$P_VAL_adj <- p.adjust(Indel_G$P_VAL, method = "bonferroni")
}

############
#### Merge indels and SNPs
############

#### TO BE DONE
VARIANTS <- rbind(SNP,Indel_A,Indel_G)

## Fisher strand filter
VARIANTS$FISHER <- apply(VARIANTS,1,function(x) fisher.test(matrix(c(as.numeric(x["Ref_fwd"]), as.numeric(x["Alt_fwd"]), as.numeric(x["Ref_rev"]), as.numeric(x["Alt_rev"])), nrow = 2))[[1]])
VARIANTS$FISHER <- p.adjust(VARIANTS$FISHER, method = "bonferroni")


## DP1-DP2 strand filter
#VARIANTS$FISHER_RATIO <- apply(VARIANTS,1,function(x) fisher.test(matrix(c(as.numeric(x["DP_DP1"]), as.numeric(x["DP_DP2"]), as.numeric(x["ALT_COUNT_DP1"]), as.numeric(x["ALT_COUNT_DP2"])), nrow = 2))[[1]])
#VARIANTS$P_RATIO <- pbinom(VARIANTS$ALT_COUNT_DP2,VARIANTS$ALT_COUNT,prob = MU, lower.tail = T)

# Distance to somatic like mutations (those ones that are outside of the Error Distribution)
SOMATIC_LIKE <- VARIANTS[VARIANTS$P_VAL_adj < 0.1,]

VARIANTS <- dist_SNP(VARIANTS,SOMATIC_LIKE)

CONTROL <- NULL


############
#### Merging results with normal tissue sites
############

# Loading normal file 

if (!is.null(ND)){
  ND <- args$normal_file
  
  ND <- as.data.frame(fread(ND))
  
  ND$AB <- ND$ALT_COUNT/ND$DP_HQ

  # Getting fastly alpha and beta for beta-binamial distribution from normal tissue
  FIT_DP1.ND <- MODEL_BETABIN_DP_HQ(ND[ND$AB < 0.1,])

  a_DP1.ND <- FIT_DP1.ND[[1]]
  b_DP1.ND <- FIT_DP1.ND[[2]]

  #### FOR TC group
  ND$P_VAL <- P_VAL(ND, a_DP1.ND, b_DP1.ND)
  
  ND_NO_VARIANTS <- ND[ND$ALT_COUNT == 0,]
  ND_VARIANTS <- ND[ND$ALT_COUNT > 0,]
  
  ND_NO_VARIANTS$nP_VAL_adj <- 1
  ND_VARIANTS$nP_VAL_adj <- p.adjust(ND_VARIANTS$P_VAL, method = "bonferroni")
  
  ND <- rbind(ND_NO_VARIANTS,ND_VARIANTS)
  
  # Keep only columns that we are interested
  nd <- ND[,c("CHROM","POS", "DP_HQ","ALT_COUNT", "nP_VAL_adj")]
  colnames(nd) <- c("CHROM","POS","DP_HQ_N","ALT_COUNT_N","nP_VAL_adj")
  
  nd$AB_N <- nd$ALT_COUNT_N / nd$DP_HQ_N
  
  # Merge with posible somatic variants (VARIANTS)
  nd_temp <- nd[nd$POS %in% VARIANTS$POS,]
  
  FINAL <- merge(VARIANTS,nd_temp, by=c("CHROM","POS"),all.x = T)
} else {
  FINAL <- VARIANTS
}

############
#### Printing results
############

write.table(x=FINAL, file=OUTF, sep='\t', quote=F, row.names=F,col.names=T)
