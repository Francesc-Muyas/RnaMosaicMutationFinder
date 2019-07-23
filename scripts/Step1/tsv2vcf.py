#!/usr/bin/env python
import argparse
import time

parser = argparse.ArgumentParser(description='Getting barcodes in fastq file and labels')
parser.add_argument('-i', '--infile', type=str, help='Tsv table', required= True)
parser.add_argument('-tID', '--tumorid', type=str, default='Tumor', help='Tumor sample id', required= False)
parser.add_argument('-nID', '--normalid', type=str, default='Normal', help='Normal sample id', required= False)
parser.add_argument('-ref', '--reference', type=str, help='Reference fastq file which table was build', required= True)
parser.add_argument('-o', '--outfile', type=str, help='Vcf output file', required= True)
parser.add_argument('-cov', '--min_COV', type=int, default=10, help='Minimum Coverage', required= False)
parser.add_argument('-ac', '--min_AC', type=int, default=3, help='Minimum reads supporting alternative allele', required= False)
parser.add_argument('-variant_dist', '--min_DIST', type=int, default=20, help='Minimum distance allowed between variants (to avoid clustered errors)', required= False)
parser.add_argument('-dup1', '--duplicate1_FILTER', type=float, default=0, help='Minimum (mean) number of duplicates 1', required= False)
parser.add_argument('-dup2', '--duplicate2_FILTER', type=float, default=0, help='Minimum (mean) number of duplicates 2', required= False)
parser.add_argument('-str', '--strand', type=int, choices = [0,1], default=1, help='Strand bias test (Fisher test). 0 for turn it off', required= False)
parser.add_argument('-af', '--min_AF', type=float, default=0, help='Minimum allele frequency allowed', required= False)
parser.add_argument('-end', '--end_read_filter', choices = [0,1],  type=int, default = 1 , help='Filtering for variants found at end or begining of the read. 0 for turn it off. Default: Activated', required= False)
parser.add_argument('-som', '--somatic_type', choices = ['single','paired'],  type=str, default = 'paired' , help='If analysis is based on two paired samples (paired) or based on only one sample (single). Default: paired', required= False)

args = parser.parse_args()


FILE = args.infile
SAMPLE = args.tumorid
NORMAL = args.normalid
REF = args.reference
out = args.outfile

min_COV = args.min_COV
min_AC = args.min_AC

end = args.end_read_filter
min_DIST = args.min_DIST
DUP1_FILTER = args.duplicate1_FILTER
DUP2_FILTER = args.duplicate2_FILTER
strand = args.strand
AF = args.min_AF

OUT_vcf = open(out,'w')


###### GETTING HEADER
date = time.strftime("%d/%m/%Y")## dd/mm/yyyy format

VCF_format="##fileformat=VCFv4.1"
DATE="##fileDate=%s" % date
source="##source=CRG_UKT_somatic_variant_calling"
reference="##reference=%s" % REF
CONCEPTS="""##INFO=<ID=Snp_Dist,Number=1,Type=Integer,Description="Distance to the closest SNP">
##FILTER=<ID=PASS,Description="Passed filter">
##FILTER=<ID=Germline,Description="Germline variant">
##FILTER=<ID=No_germline_info,Description="Not germline information">
##FILTER=<ID=Low_COV,Description="Low coverage">
##FILTER=<ID=Strand_imbalanced,Description="All alternative reads found in only one strand">
##FILTER=<ID=Low_AC,Description="Less than defined minimum of alternative counts">
##FILTER=<ID=Clustered_Variant,Description="Clustered variants">
##FILTER=<ID=Error,Description="Alternative counts inside the expected error rate distribution">
##FILTER=<ID=LQ_call,Description="Low confidence call. Very close to the error rate distribution">
##FILTER=<ID=DUP1_filter,Description="Alternative counts not supported with enough duplicates (based on barcode 1)">
##FILTER=<ID=DUP2_filter,Description="Alternative counts not supported with enough duplicates (based on barcode 2)">
##FILTER=<ID=ABB_biased,Description="Prone to systematic error based on ABB score">
##FILTER=<ID=Fisher_Strand,Description="Strand bias based on fisher test">
##FILTER=<ID=End_read,Description="More than 50 % of reads found at last of first read of the read">
##FILTER=<ID=Germline_contamination,Description="Reads supporting alternative alleles in germline sample">
##FILTER=<ID=Low_qual_pos,Description="Position enriched with too many low quality bases">
##FILTER=<ID=LQ_alt_counts,Description="Position enriched with too many low quality alternative bases">
##FILTER=<ID=Variant_contamination,Description="Reads supporting third alleles in tumor sample">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed"
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=AB,Number=1,Type=Float,Description="Allele balance">
##FORMAT=<ID=QUAL,Number=1,Type=Float,Description="Sum of alterantive qualities">
##FORMAT=<ID=DUP1,Number=1,Type=Float,Description="Mean of duplicates of all alternative counts (based on barcode 1)">
##FORMAT=<ID=DUP1,Number=1,Type=Float,Description="Mean of duplicates of all alternative counts (based on barcode 2)">
##FORMAT=<ID=Strand,Number=2,Type=String,Description="Alleles in strands, Reference forward, Reference reverse, Alternative forward, Alternative reverse">
##FORMAT=<ID=FS,Number=1,Type=Float,Description="Fisher strand test">
##FORMAT=<ID=VCB,Number=1,Type=Integer,Description="Variant Count bias, number of other different alternative alleles found">
##FORMAT=<ID=Perror,Number=1,Type=Float,Description="P-value after correction to belong to the error rate distribution (beta binomial distribution estimated in ND sample) - After bonferroni correction">"""
INFILE="##Input file:%s" % FILE
sample_name="##Tumor sample:%s" % SAMPLE
Parameters="##Parameters for filtering = min_COV: %s; min_AC: %s; min_DIST: %s; DUP1_FILTER: %s; DUP2_FILTER: %s; strand: %s; min_AF: %s; end_read_filter: %s; somatic_type: %s"  % (min_COV, min_AC, min_DIST, DUP1_FILTER, DUP2_FILTER, strand, AF, end, args.somatic_type)




## "##" VCF header
OUT_vcf.write(VCF_format + '\n')
OUT_vcf.write(DATE + '\n')
OUT_vcf.write(source + '\n')
OUT_vcf.write(reference + '\n')
OUT_vcf.write(CONCEPTS + '\n')
#OUT_vcf.write(INFILE + '\n')
#OUT_vcf.write(sample_name + '\n')
OUT_vcf.write(Parameters + '\n')

## "#" VCF row
if ( args.somatic_type == 'paired'):
    VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', SAMPLE, NORMAL]
else:
    VCF_HEADER = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', SAMPLE]

with open(FILE) as f1:
    for i in f1:
        line  =  i.rstrip('\n')
        
        if line.startswith('CHROM'):
            header_file = line
            form = header_file.split('\t')
            
            if ("CHROM" in form):
                CHROMi=[i for i, x in enumerate(form) if x == "CHROM"][0]
            else:
                print "CHROM not present in vcf"
                break

            if ("POS" in form):
                POSi=[i for i, x in enumerate(form) if x == "POS"][0]
            else:
                print "POS not present in vcf"
                break

            if ("REF" in form):
                REFi=[i for i, x in enumerate(form) if x == "REF"][0]
            else:
                print "REF not present in vcf"
                break

            if ("ALT" in form):
                ALTi=[i for i, x in enumerate(form) if x == "ALT"][0]
            else:
                print "ALT not present in vcf"
                break

            if ("Alt_qual" in form):
                QUALi=[i for i, x in enumerate(form) if x == "Alt_qual"][0]
            else:
                print "Alt_qual not present in vcf"
                break

            if ("DIST" in form):
                DISTi=[i for i, x in enumerate(form) if x == "DIST"][0]
            else:
                print "DIST not present in vcf"
                break
            
            if ("Ref_fwd" in form):
                Ref_fwdi=[i for i, x in enumerate(form) if x == "Ref_fwd"][0]
            else:
                print "Ref_fwd not present in vcf"
                break
 
            if ("Ref_rev" in form):
                Ref_revi=[i for i, x in enumerate(form) if x == "Ref_rev"][0]
            else:
                print "Ref_rev not present in vcf"
                break

            if ("Alt_fwd" in form):
                Alt_fwdi=[i for i, x in enumerate(form) if x == "Alt_fwd"][0]
            else:
                print "Alt_fwd not present in vcf"
                break

            if ("Alt_rev" in form):
                Alt_revi=[i for i, x in enumerate(form) if x == "Alt_rev"][0]
            else:
                print "Alt_rev not present in vcf"
                break
            
            if ("Start_End_count" in form):
                Start_end_counti=[i for i, x in enumerate(form) if x == "Start_End_count"][0]
            else:
                print "Start_End_count not present in vcf"
                break
            
            if ("FISHER" in form):
                FISHERi=[i for i, x in enumerate(form) if x == "FISHER"][0]
            else:
                print "FISHER not present in vcf"
                break

            if ("Lq_alt_count" in form):
                Lq_alt_counti=[i for i, x in enumerate(form) if x == "Lq_alt_count"][0]
            else:
                print "Lq_alt_count not present in vcf"
                break

            if ("Digit_1" in form):
                Digit1i=[i for i, x in enumerate(form) if x == "Digit_1"][0]
            else:
                print "Digit_1 not present in vcf"
                break

            if ("Digit_2" in form):
                Digit2i=[i for i, x in enumerate(form) if x == "Digit_2"][0]
            else:
                print "Digit_2 not present in vcf"
                break
            
            if ("VariantCountBias" in form):
                VariantCountBiasi=[i for i, x in enumerate(form) if x == "VariantCountBias"][0]
            else:
                print "VariantCountBias not present in vcf"
                break            
            
            if ("AB" in form):
                ABi=[i for i, x in enumerate(form) if x == "AB"][0]
            else:
                print "AB not present in vcf"
                break
            
            # Alternative counts
            if ("ALT_COUNT" in form):
                ALT_COUNTi=[i for i, x in enumerate(form) if x == "ALT_COUNT"][0]
            else:
                print "ALT_COUNT not present in vcf"
                break
            
            # Coverage 
            if ("DP" in form):
                DPi=[i for i, x in enumerate(form) if x == "DP"][0]
            else:
                print "DP not present in vcf"
                break
            
            if ("DP_HQ" in form):
                DP_HQi=[i for i, x in enumerate(form) if x == "DP_HQ"][0]
            else:
                print "DP_HQ not present in vcf"
                break

            if ("P_VAL_adj" in form):
                P_VAL_adji=[i for i, x in enumerate(form) if x == "P_VAL_adj"][0]
            else:
                print "P_VAL_adj not present in vcf"
                break
            
            if ( args.somatic_type == 'paired'):
                ## Paired
                # Allele balance
                
                if ("AB_N" in form):
                    AB_NDi=[i for i, x in enumerate(form) if x == "AB_N"][0]
                else:
                    print "AB_N not present in vcf"
                    break
    
                if ("ALT_COUNT_N" in form):
                    ALT_COUNT_Ni=[i for i, x in enumerate(form) if x == "ALT_COUNT_N"][0]
                else:
                    print "ALT_COUNT_N not present in vcf"
                    break
                
                if ("DP_HQ_N" in form):
                    DP_HQ_Ni=[i for i, x in enumerate(form) if x == "DP_HQ_N"][0]
                else:
                    print "DP_HQ_N not present in vcf"
                    break
                
                if ("nP_VAL_adj" in form):
                    nP_VAL_adji=[i for i, x in enumerate(form) if x == "nP_VAL_adj"][0]
                else:
                    print "nP_VAL_adj not present in vcf"
                    break

            VCF_HEADER_l = "#"+'\t'.join(VCF_HEADER)            
            
            OUT_vcf.write(VCF_HEADER_l + '\n')
            
        else:
            info = line.split("\t")
            ## COMMON VCF columns (FIRST COLUMNS)
            CHROM = info[CHROMi]
            POS = info[POSi]
            ID = '.'
            REF = info[REFi]
            ALT = info[ALTi]
            QUAL = '.'
                        
                            
            COMMON = [CHROM, POS, ID, REF, ALT, QUAL]
            
            ## INFO fields
            DIST = info[DISTi]

            INFO = ["Snp_Dist="+str(DIST)]
            
            
            ## FORMAT
            FORMAT = ["GT","AD","DP", "AB", "QUAL", "DUP1", "DUP2", "Strand", "FS", "VCB", "Perror"]
            
            if ( args.somatic_type == 'paired'):
                ## Normal 
                AB_ND = info[AB_NDi]
                ALT_COUNT_N = info[ALT_COUNT_Ni]
                DP_HQ_N = info[DP_HQ_Ni]
                
                # P-value corrected that show significance in Error distribution
                Pval_ND = info[nP_VAL_adji]
                if (Pval_ND != "NA"):
                    Pval_ND = float(round(float(Pval_ND),3))
                else:
                    Pval_ND = '.'
                
                # Checking if we have Germline info
                if (AB_ND != "NA"):
                    AB_ND = float(round(float(AB_ND),3))
                else:
                    AB_ND = '.'
                
                # Fast genotyping of Germline
                if (AB_ND != "."):
                    if (AB_ND > 0.15):
                        GT_ND = "0/1"
                    elif (AB_ND > 0.75):
                        GT_ND = "1/1"
                    else:
                        GT_ND = "0/0"
                else:
                    GT_ND = "./."
                
                # AD for normal
                if (AB_ND != "."):
                    AD_ND = str(int(DP_HQ_N) - int(ALT_COUNT_N)) + ',' + str(ALT_COUNT_N)
                else:
                    AD_ND = ".,."
                
                sample_info_N = [GT_ND, str(AD_ND), str(DP_HQ_N), str(AB_ND), ".", ".", ".", ".", ".", ".", str(Pval_ND)]
            
            ## Tumor
            QUAL = info[QUALi]
            VariantCountBias = info[VariantCountBiasi]
            ALT_COUNT = int(info[ALT_COUNTi])
            DP_HQ = info[DP_HQi]
            DP = info[DPi]
            AB = round(float(info[ABi]),3)
            Ref_fwd = int(info[Ref_fwdi])
            Ref_rev = int(info[Ref_revi])
            REF_C = int(Ref_fwd) + int(Ref_rev)
            Alt_fwd = int(info[Alt_fwdi])
            Alt_rev = int(info[Alt_revi])
            Digit1 = info[Digit1i]
            Digit2 = info[Digit2i]
            Lq_count = int(info[Lq_alt_counti])
            
            # P-value corrected that show significance in Error distribution
            Pval_TD = info[P_VAL_adji]
            if (Pval_TD != "NA"):
                Pval_TD = float(round(float(Pval_TD),3))
            else:
                Pval_TD = '.'
            
            
            SEC = info[Start_end_counti]
            if (SEC == '.'):
                SEC = 0
            Start_end_count = int(SEC)
            
            Start_end_ratio = float(Start_end_count)/ALT_COUNT
            FISHER = round(float(info[FISHERi]),4)

            Ratio_fwd_alt = round(float(Alt_fwd)/ALT_COUNT,4)
            Ratio_rev_alt = round(float(Alt_rev)/ALT_COUNT,4)
            Min_Ratio = min(Ratio_fwd_alt,Ratio_rev_alt)
            AD = str(REF_C) + ',' + str(ALT_COUNT)
             
            # Genotyping Tumor
            if (float(Pval_TD) < 0.1 and float(AB) >= 0.75):
                GT = "1/1"
            elif (float(Pval_TD) < 0.1 and float(AB) < 0.75):
                GT = "0/1"
            else:
                GT = "0/0"
                        
            sample_info = [GT, str(AD), str(DP_HQ), str(AB), str(QUAL), str(Digit1), str(Digit2), str(Ref_fwd) + "-" + str(Ref_rev) + "-" + str(Alt_fwd) + "-" + str(Alt_rev), str(FISHER), str(VariantCountBias), str(Pval_TD)]
            
            if (float(AB) > 0 ):
                filter_criteria = []
                
                if ( args.somatic_type == 'paired'):
                    if (AB_ND != "." and GT_ND != "0/0"):
                        filter_criteria.append("Germline")
                    
                    if (AB_ND != "." and ((int(ALT_COUNT) < 2*int(ALT_COUNT_N)) or Pval_ND < 0.1)):
                        filter_criteria.append("Germline_contamination")
                                   
                    if (AB_ND == "."):
                        filter_criteria.append("No_germline_info")

                if (float(Pval_TD) > 0.1):
                    filter_criteria.append("Error")
                                        
                if (float(Start_end_ratio) > 0.5 and end == 1):
                    filter_criteria.append("End_read")

                if (float(Lq_count)/(Lq_count + ALT_COUNT) > 0.25):
                    filter_criteria.append("LQ_alt_counts")
                
                if (float(AB) < AF):
                    filter_criteria.append("Low_AF")
                    
                if (min(int(Alt_fwd), int(Alt_rev)) == 0 and min(int(Ref_fwd),int(Ref_rev)) != 0 and strand == 1):
                    filter_criteria.append("Strand_imbalanced")
                
                if (int(DP_HQ) < min_COV):
                    filter_criteria.append("Low_Cov")
                    
                if (int(ALT_COUNT) < int(min_AC)):
                    filter_criteria.append("Low_AC")
                    
                if (DIST != "Inf" and int(DIST) < min_DIST):
                    filter_criteria.append("Clustered_Variant")

                if (int(DP_HQ)/float(DP) < 0.80):
                    filter_criteria.append("Low_qual_pos")

                if (float(Digit1) < DUP1_FILTER):
                    print Digit1, DUP1_FILTER
                    filter_criteria.append("DUP1_filter")

                if (float(Digit2) < DUP2_FILTER):
                    print Digit2, DUP2_FILTER
                    filter_criteria.append("DUP2_filter")
                                  
                if (int(VariantCountBias)/float(ALT_COUNT) > 0.1):
                    filter_criteria.append("Variant_contamination")
                    
                #if (float(FISHER) < 0.001 and strand == 1):
                #    filter_criteria.append("Fisher_Strand")
                              
                if (len(filter_criteria) == 0):
                    FILTER = "PASS"
                else:
                    FILTER = ';'.join(filter_criteria)
                
                if ( args.somatic_type == 'paired'):
                    FINAL_LINE=['\t'.join(COMMON), FILTER, ';'.join(INFO), ':'.join(FORMAT), ':'.join(sample_info), ':'.join(sample_info_N)]
                else:
                    FINAL_LINE=['\t'.join(COMMON), FILTER, ';'.join(INFO), ':'.join(FORMAT), ':'.join(sample_info)]
            
                OUT_vcf.write('\t'.join(FINAL_LINE) + '\n')
            
OUT_vcf.close()
            
            
            
            