#!/usr/bin/python

import argparse
import numpy

def coverage(quality_sequence, minBQ):
    COV = len(quality_sequence)
    HQ = len([x for x in range(0,len(quality_sequence)) if (ord(quality_sequence[x])-33) >= minBQ])
    DP2 = len([x for x in range(0,len(quality_sequence)) if (ord(quality_sequence[x])-33) >= 20])

    return COV, HQ, DP2

def alt_count_2dp(digit):
    DP2 = len([x for x in range(0,len(digit)) if digit[x] >= 2])
    return DP2

def similarity(a, b):
    ratio = SequenceMatcher(None, a, b).ratio()
    return ratio

def alt_qual_digit(variant, quality, DIG1,DIG2):
    if variant.upper() in DIG1:
	DIG1[variant.upper()].append(int(quality/10))
    	DIG2[variant.upper()].append(int(quality%10))
    else:
	DIG1[variant.upper()] = list()
	DIG2[variant.upper()] = list()
	
	DIG1[variant.upper()].append(quality/10)
    	DIG2[variant.upper()].append(quality%10)
    
    return (DIG1, DIG2)

def FWD_RVS_STRAND(a):
    STRAND = "forward" if any(map(str.isupper, a)) else "reverse"
    return STRAND


def ReverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in reversed(seq)])
	
def Complement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in seq])


def keywithmaxval(d):
     v=list(d.values())
     k=list(d.keys())
     return k[v.index(max(v))]
    

def pileup_INFO(line, indel_countdown):
	line = line.rstrip('\n')
	
	# The are some lines in some pileup files that do not have the expected number of fields
	if len(line.split("\t")) != 6:
	    return ".", indel_countdown-1
	
	locus, pos, ref_base, cov, bases, qualities = line.split("\t")
	
	keep_REF = ''
	keep_base_quality = ''
	
	pos = int(pos)
	cov = int(cov)
		
	b_pos = 0
	q_pos = 0
	mm_count = 0
	mm_count_lq = 0
	ref_count = 0
	ref_count_fwd = 0
	ref_count_rev = 0
	ref_count_lq = 0
	N_count = 0
	alt_counts = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0 }
	alt_quals  = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0 }
	alt_quals_digit1 = {}
	alt_quals_digit2 = {}
	
	del_quals = {}
	del_quals_digit1 = {}
	del_quals_digit2 = {}
	
	ins_quals = {}
	ins_quals_digit1 = {}
	ins_quals_digit2 = {}
	
	alt_ends = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0 }
	alt_starts = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0 }
	alt_overlapped = { 'A': 0, 'C': 0, 'G': 0, 'T': 0, 'a': 0, 'c': 0, 'g': 0, 't': 0 }
	indel_counts = {'del': 0, 'ins': 0}
	
	insertion=dict()
	insertion_FWD=dict()
	insertion_RVS=dict()

	deletion=dict()
	deletion_FWD=dict()
	deletion_RVS=dict()
	
	#total_pos += 1
	END = 0
	START= 0

	### Parse pileup line
	while b_pos < len(bases):
		STRAND = 0
		base = bases[b_pos]
		    
		# Match
		if base in '.,':
			base_quality = ord(qualities[q_pos]) - 33
			keep_REF = base
			keep_base_quality = base_quality
			if base_quality >= minBQ:
			    ref_count += 1
			    
			    if base == ".":	
				ref_count_fwd += 1
				
			    elif base == ",":
				ref_count_rev += 1
			    
			else:
				ref_count_lq += 1
				#total_lq_rf += 1
			
			b_pos += 1
			q_pos += 1
			
			END = 0
			START = 0
			
		# Mismatch
		elif base in 'ATCGatcg':
			base_quality = ord(qualities[q_pos]) - 33
			if base_quality >= minBQ:
				alt_counts[base] += 1
				alt_quals[base] += base_quality
				alt_quals_digit1, alt_quals_digit2 = alt_qual_digit(base, base_quality, alt_quals_digit1, alt_quals_digit2)
				mm_count += 1
					
				if END > 0:
				    alt_ends[base] += 1
				    
				elif START > 0:
				    alt_starts[base] += 1
				
			else:
				mm_count_lq += 1
				#total_lq_mm += 1
			
			b_pos += 1
			q_pos += 1
			    
			END = 0
			START = 0
			
		# Ambiguous base (N)
		elif base in 'Nn':
			keep_REF = base
			b_pos += 1
			q_pos += 1
			N_count += 1
			
			END = 0
			START = 0
			
		# End of read sign
		elif base == '$':
			b_pos += 1
			END = 1
			START = 0
			
		# Start of read sign, followed by MapQ ascii value
		elif base == '^':
			b_pos += 2
			START = 1
			END = 0
			
		# Deletion
		elif base == '-':
			### Before indels, there is always a ref base where the indel is happening. This has not to be count as reference count because is refering to the indel.
			
			if (keep_REF == "." and keep_base_quality >= minBQ):
			    ref_count = ref_count -1 
			    ref_count_fwd = ref_count_fwd - 1
			elif (keep_REF == "," and keep_base_quality >= minBQ):
			    ref_count = ref_count -1 
			    ref_count_rev = ref_count_rev - 1
			else:
			    # get deletion length and the deletion
			    b_pos += 1
			    indel_len = bases[b_pos]
			    
			    while bases[b_pos+1].isdigit():
				indel_len=indel_len+bases[b_pos+1] ##We concatenate the numbers to create a number of more digits 
				b_pos += 1
			    
			    b_pos += 1
			    
			    
			    # jump over deletion
		    	    b_pos += int(indel_len)
			    
			    END = 0
			    START = 0
			
			    continue
			
			
			indel_counts['del'] += 1
			
			# get deletion length and the deletion
			b_pos += 1
			indel_len = bases[b_pos]
			
			
			while bases[b_pos+1].isdigit():
				indel_len=indel_len+bases[b_pos+1] ##We concatenate the numbers to create a number of more digits 
				b_pos += 1
			
			indel_len=int(indel_len)
			#while bases[b_pos].isdigit():
			#	indel_len += bases[b_pos]
			#	b_pos += 1
			
			del_observed= bases[(int(b_pos)+1):(int(b_pos)+int(indel_len)+1)]
			alt_quals_digit1, alt_quals_digit2 = alt_qual_digit(del_observed, base_quality, alt_quals_digit1, alt_quals_digit2)
			
			
			STRAND = FWD_RVS_STRAND(del_observed) ## To know the strand (Foward or reverse)
			
			del_observed=del_observed.upper() ## With this command we cannot ditinguish between strands
			
			del_quals_digit1, del_quals_digit2 = alt_qual_digit(del_observed, keep_base_quality, del_quals_digit1, del_quals_digit2)
			
			if del_observed in deletion:
				deletion[del_observed] += 1
				del_quals[del_observed] += keep_base_quality
				
				if STRAND == "forward":
					if del_observed in deletion_FWD:
						deletion_FWD[del_observed] += 1
					else:
						deletion_FWD[del_observed] = 1
				
				elif STRAND == "reverse":
					if del_observed in deletion_RVS:
						deletion_RVS[del_observed] += 1
					else:
						deletion_RVS[del_observed] = 1
			else:
				deletion[del_observed] = 1
				del_quals[del_observed] = keep_base_quality
				
				if STRAND == "forward":
					deletion_FWD[del_observed] = 1
					deletion_RVS[del_observed] = 0
				
				elif STRAND == "reverse":
					deletion_RVS[del_observed] = 1
					deletion_FWD[del_observed] = 0
					
			b_pos += 1
			
			
			# jump over deletion
			b_pos += int(indel_len)
			
			END = 0
			START = 0


		# Insertion
		elif base == '+':
			### Before indels, there is always a ref base where the indel is happening. This has not to be count as reference count because is refering to the indel.
			 
			if (keep_REF == "." and keep_base_quality >= minBQ):
			    ref_count_fwd = ref_count_fwd - 1
			    ref_count = ref_count -1
			elif (keep_REF == "," and keep_base_quality >= minBQ):
			    ref_count_rev = ref_count_rev - 1
			    ref_count = ref_count -1
			else:
			    b_pos += 1
			    indel_len = bases[b_pos]
		    	    
			    while bases[b_pos+1].isdigit():
				    indel_len=indel_len+bases[b_pos+1] ##We concatenate the numbers to create a number of more digits 
				    b_pos += 1
				    
			    indel_len=int(indel_len)
			    
			    b_pos += 1
			    
			    b_pos += int(indel_len)
			        
			    END = 0
			    START = 0
			    
			    continue
			
			indel_counts['ins'] += 1
			
			# get insertion length
			b_pos += 1
			indel_len = bases[b_pos]
			
			while bases[b_pos+1].isdigit():
				indel_len=indel_len+bases[b_pos+1] ##We concatenate the numbers to create a number of more digits 
				b_pos += 1
			
			indel_len=int(indel_len)
			
			ins_observed= bases[(int(b_pos)+1):(int(b_pos)+int(indel_len)+1)]
			
			ins_quals_digit1, ins_quals_digit2 = alt_qual_digit(ins_observed, keep_base_quality, ins_quals_digit1, ins_quals_digit2)

			STRAND = FWD_RVS_STRAND(ins_observed)
			
			ins_observed=ins_observed.upper()
			
			if ins_observed in insertion:
				insertion[ins_observed] += 1
				ins_quals[ins_observed] += keep_base_quality
				
				if STRAND == "forward":
					if ins_observed in insertion_FWD:
						insertion_FWD[ins_observed] += 1
					else:
						insertion_FWD[ins_observed] = 1
				
				elif STRAND == "reverse":
					if ins_observed in insertion_RVS:
						insertion_RVS[ins_observed] += 1
					else:
						insertion_RVS[ins_observed] = 1
			else:
				insertion[ins_observed] = 1
				ins_quals[ins_observed] = keep_base_quality

				if STRAND == "forward":
					insertion_FWD[ins_observed] = 1
					insertion_RVS[ins_observed] = 0
				elif STRAND == "reverse":
					insertion_RVS[ins_observed] = 1
					insertion_FWD[ins_observed] = 0
			

			b_pos += 1
			
			#while bases[b_pos].isdigit():
			#	indel_len += bases[b_pos]
			#	b_pos += 1

			# jump over insertion
			b_pos += int(indel_len)

			END = 0
			START = 0

		# Deletion placeholder
		elif base == '*':
			b_pos += 1

		# Crap
		else:
			print "ERROR:", locus, pos, b_pos, base
			exit(1)


	### mask 10bp behind indels
	if(indel_counts['del'] + indel_counts['ins'] > 5):
		indel_countdown = 0
	indel_countdown += 1


	#################
	# NEW IDEA. "." or "," before the indels explain quality of the indels, they are not ref base counts.
	#################
	
	#ref_count = ref_count - indel_counts['del'] - indel_counts['ins']
	#print ref_count, ref_count_fwd, ref_count_rev
	DP, DP_HQ, DP_DP2 = coverage(qualities, minBQ)
	
	### SNP/indel calling
	
	# sum up counts from fwd and rev strand
	alt_count_total = {}
	alt_count_total['A'] = alt_counts['A'] + alt_counts['a']
	alt_count_total['C'] = alt_counts['C'] + alt_counts['c']
	alt_count_total['G'] = alt_counts['G'] + alt_counts['g']
	alt_count_total['T'] = alt_counts['T'] + alt_counts['t']
	
	# sum up base qualities from fwd and rev strand
	alt_qual_total = {}
	alt_qual_total['A'] = alt_quals['A'] + alt_quals['a']
	alt_qual_total['C'] = alt_quals['C'] + alt_quals['c']
	alt_qual_total['G'] = alt_quals['G'] + alt_quals['g']
	alt_qual_total['T'] = alt_quals['T'] + alt_quals['t']
	
	alt_call  = keywithmaxval(alt_count_total)
	alt_count = alt_count_total[alt_call]
	
	# if there are alternative allele counts, we have to decide between snps and indels
	if ((sum(alt_count_total.values()) + sum(indel_counts.values())) > 1):
	    
	    # We have to decide if to call SNP or Indel looking at the counts
	    # Get counts for most common SNP, Indel and Deletion
	    VARIANT_TYPE = {}
	    
	    # Snp
	    SNP_count = alt_count
	    VARIANT_TYPE['SNP'] = SNP_count
	    
	    # Delition
	    if (indel_counts['del'] > 0):
		del_call = keywithmaxval(deletion)
		Del_count = deletion[del_call]
	    else:
		Del_count = 0
	    VARIANT_TYPE['Del'] = Del_count
    
	    # Insertion
	    if (indel_counts['ins'] > 0):
		ins_call = keywithmaxval(insertion)
		Ins_count = insertion[ins_call]
	    else:
		Ins_count = 0
	    VARIANT_TYPE['Ins'] = Ins_count
	    
	    # We take the max key
	    MAX = keywithmaxval(VARIANT_TYPE)
	    if ( MAX == 'SNP'):
		types = 'SNP'
	    else:
		# If SNP counts are equal to indel count. We take SNP as default
		if (VARIANT_TYPE[MAX] > VARIANT_TYPE['SNP']):
		    types = 'INDEL'
		else:
		    types = 'SNP'
	else:
	    types = 'SNP'
	
	if types == "SNP":

	
		# Identify major alternative allele
		# alt_call =  max(d, key=d.get)
		alt_call  = keywithmaxval(alt_count_total)
		alt_count = alt_count_total[alt_call]
		alt_qual  = alt_qual_total[alt_call]
		alt_count_fwd = alt_counts[alt_call.upper()]
		alt_count_rev = alt_counts[alt_call.lower()]
		alt_count_end = alt_ends[alt_call.upper()] + alt_ends[alt_call.lower()] 
		alt_count_start = alt_starts[alt_call.upper()] + alt_starts[alt_call.lower()] 
		alt_count_overlapped = alt_overlapped[alt_call.upper()] + alt_overlapped[alt_call.lower()]
		

		OTHER_READS = (mm_count + indel_counts['ins'] + indel_counts['del']) - alt_count
		
		# Call SNP using simple thresholds
		if ((alt_count >=1) and ((alt_count + ref_count)>= min_DP) and (alt_qual/alt_count >= minBQ) ): #and indel_countdown <= 9):
			
			ALT_DP2 = alt_count_2dp(alt_quals_digit1[alt_call.upper()])
			
			DIGIT_1 = round(numpy.mean(alt_quals_digit1[alt_call.upper()]),2)
			DIGIT_2 = round(numpy.mean(alt_quals_digit2[alt_call.upper()]),2)

			line = (locus, str(pos), str(ref_base.upper()), alt_call, "SNP", str(DP), str(DP_HQ), str(DP_DP2), str(alt_count), str(ALT_DP2), str(alt_qual), str(DIGIT_1), str(DIGIT_2), str(ref_count_fwd), str(ref_count_rev), str(alt_count_fwd), str(alt_count_rev), str(alt_count_end + alt_count_start), str(OTHER_READS), str(N_count), str(mm_count_lq))

                        line = "\t".join(line)

			return line,indel_countdown		

		elif ((alt_count == 0) and ((alt_count + ref_count)>= min_DP)): #and indel_countdown <= 9):
			
			ALT_DP2 = 0
			
			line = (locus, str(pos), str(ref_base.upper()), ".", ".", str(DP), str(DP_HQ), str(DP_DP2), str(alt_count), str(ALT_DP2), ".", ".", ".", str(ref_count_fwd), str(ref_count_rev), str(0), str(0), str(0), str(OTHER_READS), str(N_count), str(mm_count_lq))
                        line = "\t".join(line)


                        return line, indel_countdown
                
                else:
                    return "LOW_QUAL",indel_countdown-1
                
	elif types == "INDEL":

		if len(deletion) > 0 and len(insertion) >0 and DP_HQ >= min_DP:
			del_call = keywithmaxval(deletion)
			del_qual = del_quals[del_call]
			DIGIT_1d = round(numpy.mean(del_quals_digit1[del_call.upper()]),2)
			DIGIT_2d = round(numpy.mean(del_quals_digit2[del_call.upper()]),2)
			
			ins_call = keywithmaxval(insertion)
			ins_qual = ins_quals[ins_call]
			DIGIT_1i = round(numpy.mean(ins_quals_digit1[ins_call.upper()]),2)
			DIGIT_2i = round(numpy.mean(ins_quals_digit2[ins_call.upper()]),2)
			
			if deletion[del_call] >= insertion[ins_call]:
				
				### DELETION
				
				ALT_DP2 = alt_count_2dp(del_quals_digit1[del_call.upper()])
				
				OTHER_READS = mm_count + indel_counts['ins']
				DIGIT_1d = round(numpy.mean(del_quals_digit1[del_call.upper()]),2)
				DIGIT_2d = round(numpy.mean(del_quals_digit2[del_call.upper()]),2)
				
				#### Getting ref and alternative left aligned 
				REF_base = ref_base.upper() + del_call.upper()
				DEL_call = ref_base.upper()
								
				line = (locus, str(pos), str(REF_base.upper()), DEL_call, "del", str(DP), str(DP_HQ), str(DP_DP2), str(deletion[del_call]), str(ALT_DP2), str(del_qual),  str(DIGIT_1d), str(DIGIT_2d), str(ref_count_fwd) , str(ref_count_rev), str(deletion_FWD[del_call]), str(deletion_RVS[del_call]), ".", str(OTHER_READS), str(N_count), str(mm_count_lq))
				
				line = "\t".join(line)

                                return line, indel_countdown
                            
                        else:
				### INSERTION
				
				ALT_DP2 = alt_count_2dp(ins_quals_digit1[ins_call.upper()])
				
				OTHER_READS = mm_count + indel_counts['del']
				DIGIT_1i = round(numpy.mean(ins_quals_digit1[ins_call.upper()]),2)
				DIGIT_2i = round(numpy.mean(ins_quals_digit2[ins_call.upper()]),2)

				#### Getting ref and alternative left aligned 
				REF_base = ref_base.upper()
				INS_call = ref_base.upper() + ins_call.upper()

				line = (locus, str(pos), str(REF_base.upper()), INS_call, "ins", str(DP), str(DP_HQ), str(DP_DP2), str(insertion[ins_call]), str(ALT_DP2), str(ins_qual), str(DIGIT_1i), str(DIGIT_2i), str(ref_count_fwd) , str(ref_count_rev), str(insertion_FWD[ins_call]), str(insertion_RVS[ins_call]), ".", str(OTHER_READS), str(N_count), str(mm_count_lq))
				
				line = "\t".join(line)

                                return line, indel_countdown
		
		elif len(deletion) > 0 and len(insertion) == 0 and DP_HQ >= min_DP:
			
			### DELETION
			
			ALT_DP2 = alt_count_2dp(del_quals_digit1[del_call.upper()])
			
			del_call = keywithmaxval(deletion)
			del_qual = del_quals[del_call]
			DIGIT_1d = round(numpy.mean(del_quals_digit1[del_call.upper()]),2)
			DIGIT_2d = round(numpy.mean(del_quals_digit2[del_call.upper()]),2)

			
			OTHER_READS = mm_count + indel_counts['ins']

			#### Getting ref and alternative left aligned 
			REF_base = ref_base.upper() + del_call.upper()
			DEL_call = ref_base.upper()
								
			line = (locus, str(pos), str(REF_base.upper()), DEL_call, "del", str(DP), str(DP_HQ), str(DP_DP2), str(deletion[del_call]), str(ALT_DP2), str(del_qual),  str(DIGIT_1d), str(DIGIT_2d), str(ref_count_fwd) , str(ref_count_rev), str(deletion_FWD[del_call]), str(deletion_RVS[del_call]), ".", str(OTHER_READS), str(N_count), str(mm_count_lq))

			line = "\t".join(line)

			return line, indel_countdown
							
		elif len(deletion) == 0 and len(insertion) > 0 and ((ref_count + mm_count + indel_counts['del'] + indel_counts['ins'])>= min_DP):
		    
			### INSERTION
			ALT_DP2 = alt_count_2dp(ins_quals_digit1[ins_call.upper()])
			
			ins_call = keywithmaxval(insertion)
			ins_qual = ins_quals[ins_call]
			DIGIT_1i = round(numpy.mean(ins_quals_digit1[ins_call.upper()]),2)
			DIGIT_2i = round(numpy.mean(ins_quals_digit2[ins_call.upper()]),2)

			OTHER_READS = mm_count + indel_counts['del']
			

			#### Getting ref and alternative left aligned
			REF_base = ref_base.upper()
			INS_call = ref_base.upper() + ins_call.upper()
			
			line = (locus, str(pos), str(REF_base.upper()), INS_call, "ins", str(DP), str(DP_HQ), str(DP_DP2), str(insertion[ins_call]), str(ALT_DP2), str(ins_qual), str(DIGIT_1i), str(DIGIT_2i), str(ref_count_fwd) , str(ref_count_rev), str(insertion_FWD[ins_call]), str(insertion_RVS[ins_call]), ".", str(OTHER_READS), str(N_count), str(mm_count_lq))
			
			line = "\t".join(line)

			return line, indel_countdown
                        
		elif len(deletion) == 0 and len(insertion) == 0 and ((ref_count + mm_count + indel_counts['del'] + indel_counts['ins']) >= min_DP):
			ALT_DP2 = 0
			
			OTHER_READS = mm_count + indel_counts['ins'] + indel_counts['del']

			line = (locus, str(pos), str(ref_base.upper()), ".", ".", str(DP), str(DP_HQ), str(DP_DP2), str(0), str(ALT_DP2), ".", ".", ".", str(ref_count_fwd) , str(ref_count_rev), str(0), str(0), ".", str(OTHER_READS), str(N_count), str(mm_count_lq))

			line = "\t".join(line)

                        return line, indel_countdown
                    
                else:
                    return "LOW_QUAL",indel_countdown-1
                    

#######################################################################
##	RUNNING THE CALLER
#######################################################################

### Read parameters
parser = argparse.ArgumentParser(description='Predict SNPs and indels from pileup files.')
parser.add_argument('--pileup', type=argparse.FileType('r'), required=True, dest='pileup', help='Pileup file.')
parser.add_argument('-o','--variantf', type=str, required=True, dest='variantf', help='Output file for predicted SNPs and indels.')
parser.add_argument('--minBQ', required=False, dest='minBQ', type=int, default=10, help='Minimum base quality to consider nucleotide in SNP analysis. Default = 10')
parser.add_argument('-m','--mindepth', required=False, dest='mindepth', type=int, default=1, help='Minimum minimum total depth')

### Read parameters
args = ''
try:
	args = parser.parse_args()
except IOError as io:
	print io
	sys.exit('Error reading parameters.')



pileup_file  = args.pileup
OUT_CALLING = args.variantf
minBQ = args.minBQ
min_DP = args.mindepth

indel_countdown = 1000

# Open the file to write the variant calls
OUT_CALLING = open(OUT_CALLING,'w')

### HEADER OF THE OUTPUT FILE
HEADER = ("CHROM", "POS", "REF", "ALT", "VARIANT_TYPE", "DP", "DP_HQ", "DP_DP2", "ALT_COUNT", "ALT_COUNT_DP2", "Alt_qual", "Digit_1", "Digit_2", "Ref_fwd", "Ref_rev", "Alt_fwd","Alt_rev", "Start_End_count" , "VariantCountBias", "N_count", "Lq_alt_count")
HEADER = ("\t".join(HEADER))+"\n"
OUT_CALLING.write(HEADER)

# Parse pileup file
for CHROM_POS in pileup_file:
    SPLIT = CHROM_POS.split("\t")
    
    if len(SPLIT) == 6 and int(SPLIT[3]) >= 5:

        INFO = pileup_INFO(CHROM_POS, int(indel_countdown))

	if INFO[0] != "LOW_QUAL":
            
           
	    INFO_POS = INFO[0]+"\n"
	    OUT_CALLING.write(INFO_POS)
            
        indel_countdown = INFO[1]
 
 
   
### Finish
pileup_file.close()
OUT_CALLING.close()
exit(0)