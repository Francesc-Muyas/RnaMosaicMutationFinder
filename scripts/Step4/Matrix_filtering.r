suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(argparse))


parser <- ArgumentParser()

# setting parameters
parser$add_argument("-t", "--tsv", type="character", help="Tsv file with counts", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-m", "--matrix", type="character", help="Matrix table with all tissues and individuals", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-rf", "--random_forest", type="character", help="Random forest R object", nargs=1, required=TRUE)
parser$add_argument("-b", "--black", type="character", help="Bedfile of blacklisted genes", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-e", "--exclude", type="character", help="Bed file with excluding regions (editing, systematic errors, germline sites...)", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-o", "--output_file", type="character", help="output_file prefix", metavar="file", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()
tsv <- args$tsv
Matrix <- args$matrix
random_forest <- args$random_forest
black <- args$black
exclude <- args$exclude
OUT <- args$output_file

## Primary data
# Random forest object
rf <- readRDS(random_forest)

# Bed file of blacklisted genes
BLACK <- as.data.frame(fread(black, header = F))
BLACK$V1 <- sub('chr', '', BLACK$V1)
BLACK_bed <- with(BLACK, GRanges(V1, IRanges(V2, V3)))

# Excluding bed file
EXCLUDING <- as.data.frame(fread(exclude, header = F))
EXCLUDING$V1 <- sub('chr', '', EXCLUDING$V1)
EXCLUDING_bed <- with(EXCLUDING, GRanges(V1, IRanges(V2, V3)))

##############################
## TSV random forest
##############################
TSV <- as.data.frame(fread(tsv))
colnames(TSV)[1] <- 'CHR'
TSV$CHR <- sub('chr', '', TSV$CHR)

## Columns required for random forest model
TSV <- TSV[TSV$ALT_COUNT > 3,]

# Fisher strand
TSV$FISHER_adj <- TSV$FISHER

# Alternative quality
TSV$QUAL_base <- TSV$Alt_qual/TSV$ALT_COUNT

# Strand imbalance
TSV$STRAND_REF <- rowSums(TSV[,c("Ref_fwd", "Ref_rev")] > 3)
TSV$STRAND_ALT <- rowSums(TSV[,c("Alt_fwd", "Alt_rev")] > 0)
TSV$STRAND <- (TSV$STRAND_ALT >= TSV$STRAND_REF)*1

# Black listed
TSV_bed <- with(TSV, GRanges(CHR, IRanges(POS, POS)))
TSV$BLACK <- (TSV_bed %within% BLACK_bed)*1
TSV[TSV$BLACK == 0,]$BLACK <- 2

# Random forest prediction
TSV$AB_PRED <- predict(rf, newdata = TSV, type = 'prob')[,2]

TSV$HQ_PASS <- (TSV$AB_PRED >= 0.19 &
               TSV$ALT_COUNT > 3 &
               TSV$AB >= 0.05 &
               TSV$AB < 0.75 &
               TSV$VariantCountBias < 2)*1

TSV$Soft_PASS <- TSV$FILTER 
HQ <- TSV[TSV$HQ_PASS == 1,]

##############################
## MATRIX
##############################
CALLS <- as.data.frame(fread(Matrix))
CALLS$CHR <- sub('chr', '', CALLS$CHR)

COL <- ncol(CALLS)
N_SAMPLES <- length(unique(CALLS$Id))

# Removing multiallelic sites
CALLS <- CALLS[nchar(CALLS$ALT) < 2 & CALLS$ALT != '.',]

# Creating NUM_CALL column
CALLS$NUM_CALL <- apply(CALLS, 1, function(x) sum(x[c(6:COL)] %in% c('A','C', 'T', 'G')))

# Creating NUM_FP column
CALLS$NUM_FP <- apply(CALLS, 1, function(x) sum(grepl("FP*", x[c(6:COL)], perl=TRUE)))

# Creating NUM_EXP
CALLS$NUM_EXP <- apply(CALLS, 1, function(x) sum(!is.na(x[c(6:COL)])))

# Removing multiallelic sites
CALLS <- CALLS[nchar(CALLS$ALT) < 2 & nchar(CALLS$ALT) < 2 & CALLS$ALT != '.',]

# Removing sites supported by only one FP
CALLS <- CALLS[!(CALLS$NUM_CALL == 0 & CALLS$NUM_FP < 2),]

# Editing and systematic 
CALLS_bed <- with(CALLS, GRanges(CHR, IRanges(POS, POS)))
CALLS$EXCLUDE <- (CALLS_bed %within% EXCLUDING_bed)*1

## Checking for recurrence of calls and low quality calls
TEST <- CALLS[,c("CHR","POS","Id", "NUM_CALL","NUM_FP","NUM_EXP")]
TEST$Length <- 1

# Filtering by recurrence of False Positives or low quality calls
TEST$fp <- (TEST$NUM_FP > 0)*1
FalsePositive <- aggregate(cbind(fp,Length) ~ CHR+POS, data=TEST, sum)
FalsePositiveList <- FalsePositive[FalsePositive$fp > 2,]

# Filtering for positions with recurrent 'somatic calls'
TEST$alternative <- (TEST$NUM_CALL > 0)*1
Somatic <- aggregate(cbind(alternative,Length) ~ CHR+POS, data=TEST, sum)
SomaticList <- Somatic[Somatic$alternative > 2,]

# Removing recurrent somatic and low quality calls
CALLS_filtered <- CALLS[CALLS$NUM_CALL > 0 &
                          CALLS$NUM_FP/(CALLS$NUM_CALL + CALLS$NUM_FP) <= 0.25 &
                          CALLS$NUM_CALL/CALLS$NUM_EXP < 1 &
                          CALLS$EXCLUDE != 1 &
                          !(paste(CALLS$CHR, CALLS$POS, sep='-')) %in% paste(FalsePositiveList$CHR,FalsePositiveList$POS,sep='-') &
                          !(paste(CALLS$CHR, CALLS$POS, sep='-')) %in% paste(SomaticList$CHR,SomaticList$POS,sep='-') &
                          (paste(CALLS$CHR, CALLS$POS, CALLS$Id, sep='-')) %in% paste(HQ$CHR,HQ$POS,HQ$Id,sep='-'),]



##############################
## FINAL FILES
##############################
TSV_FINAL <- TSV[paste(TSV$CHR,TSV$POS,TSV$Id,sep='-') %in% paste(CALLS$CHR, CALLS$POS, CALLS$Id, sep='-') & (TSV$FILTER == 1 | TSV$HQ_PASS == 1) ,]

# Individual calls
FINAL <- CALLS_filtered

# Ploting files
OUT1 <- OUT
OUT2 <- paste(OUT,'all_info', sep='.')
OUT3 <- paste(OUT,'all_info_filtered', sep='.')
OUT_FP <- paste(OUT,'low_quality', sep='.')
OUT_SOM <- paste(OUT,'som_recurrent', sep='.')
OUT_SINGLE <- paste(OUT,'variant_calling', sep='.')


write.table(FINAL,file=OUT1,quote=F,sep='\t',col.names=T,row.names=F)
write.table(CALLS,file=OUT2,quote=F,sep='\t',col.names=T,row.names=F)
write.table(FINAL,file=OUT3,quote=F,sep='\t',col.names=T,row.names=F)
write.table(FalsePositive,file=OUT_FP,quote=F,sep='\t',col.names=T,row.names=F)
write.table(Somatic,file=OUT_SOM,quote=F,sep='\t',col.names=T,row.names=F)
write.table(TSV_FINAL,file=OUT_SINGLE,quote=F,sep='\t',col.names=T,row.names=F)






