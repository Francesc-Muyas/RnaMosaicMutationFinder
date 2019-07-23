import pysam
import argparse


# Quality filter flag. 1 if it passes the filter and 0 if not.
def QC_read(read):
	read_out=''
	if (read.is_unmapped or read.mapq < args.minMQ or read.cigarstring == None or (read.cigarstring.count("I") + read.cigarstring.count("D")) > args.maxGAP):
		read_out = 0
	
	else:
		INDEL_SIZES=0
		if  ((read.cigarstring.count("I") + read.cigarstring.count("D")) > 0):
			
			for (cigarType,cigarLength) in read.cigar:
				try:
					if (cigarType == 1 or cigarType == 2):
						INDEL_SIZES = INDEL_SIZES + cigarLength
				except:
					continue
				
						
		NM = read.opt("NM") - INDEL_SIZES
		
		if (NM > args.maxMM):
			#print read.qname, read.cigarstring, read.opt("NM"), NM
			read_out = 0

		else:
			read_out = 1
		
	return (read_out)
	


### Read parameters
parser = argparse.ArgumentParser(description='Creates a statistic on mismatch and gap frequencies in BAM files.')
parser.add_argument('--infile', required=True, dest='infile', help='Input BAM file.')
parser.add_argument('--outfile', required=True, dest='outfile', help='Output BAM file.')
parser.add_argument('--minMQ', required=False, dest='minMQ', type=int, default=30, help='Minimum required mapping quality of aligned read. Default = 30')
parser.add_argument('--maxMM', required=False, dest='maxMM', type=int, default=4, help='Maximum number of mismatches of aligned read . Default = 4')
parser.add_argument('--maxGAP', required=False, dest='maxGAP', type=int, default=1, help='Maximum number of GAPs . Default = 1')


args = ''
try:
	args = parser.parse_args()
except IOError as io:
	print io
	sys.exit('Error reading parameters.')


### Input BAM
try:
	infile = pysam.Samfile( args.infile, "rb" )
except:
	exit("Cannot open input file.")


### Output BAM
try:
	outfile = pysam.Samfile(args.outfile, mode="wb", template = infile)
except:
	exit("Cannot open output file.")


COUNT = 0
LQ = 0
### Parse BAM file.
while 1:

	# Get read
	try:
		read = infile.next()
	except StopIteration:
		break
	
	COUNT = COUNT + 1
	
	# Filter
	READ_out = QC_read(read)
	
	if( READ_out == 1 ):
		outfile.write(read)
	else:
		#print read.mapq, read.cigarstring
		LQ = LQ + 1
		
### Finish

print "Removed {} of {} pairs due to low quality.".format(LQ, COUNT)

infile.close()
outfile.close()
exit(0)


### pysam documentation
# http://pysam.readthedocs.org/en/latest/
