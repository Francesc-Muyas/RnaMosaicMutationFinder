#!/usr/bin/env python


def read_one_line(chr_file,pos_file,cur_chr,cur_pos):
    #It checks that the current file is correct to read one more line
    #meanstaht the chr should be the same and the position should be <= to the current position
    if chr_file == cur_chr and pos_file <= cur_pos:
        return True
    else :
        return False
    
def updatefields(GT,REF, ALT, chrs,pos,i,line): #It should check that if the line is null it should just update some flag
    fields=line.split('\t')
    if line != '':

        if fields[0]==chrs[i] and  pos[i]>=int(fields[1]):
		pass
        else:
	    #print fields[1]
            chrs[i]=fields[0]
            pos[i]=int(fields[1])
	    REF[i]=fields[2]
	    ##REF = str(fields[2])
	    ALT[i]=fields[3]
            GT[i]=fields[len(fields)-1]

	    if (GT[i] == "1"):
		GT[i] = fields[3]
		
	    if(GT[i] == "0"):
		GT[i] = "."
        
	    if (GT[i] == "FP"):
		GT[i] = fields[len(fields)-1] + "." + fields[3]
	    
    return (GT,REF,ALT,chrs,pos)

def calc_low(chrs,pos,cur_chr,cur_pos): #calculate the lowest current position
    
    pos_in_cur_chr= [pos[x] for x in range(len(chrs)) if chrs[x] == cur_chr  and pos[x]>cur_pos]
	
    try:
        low_pos = min(pos_in_cur_chr)
        if low_pos <= cur_pos:
            print 'error'
            print low_pos
            raise
    except:
        print 'error'
        print chrs
        print pos
        print cur_chr
        print cur_pos
        raise
    return  low_pos


def check_chr(cur_chr,chrs):
    #it returns true if at least one chr is equal to the current chr
    for c in chrs:
        if c == cur_chr:
            return True
    return False

def update_chr(chrs):
    #retruns the current chr from the  reads of the files
    #print chrs
    #print [s for s in chrs if s]
    new_chr = [s for s in chrs if s][0]
    return new_chr

    
def write_me(chrs, pos, REF, ALT, ID, GT,cur_chr,outfile,low_pos):
    #extract the correct GT to write to the output file
    #for each of the files check those at the lowest position possible
    # and produce the output from them
    #low_pos_GT =[str(cur_chr),str(low_pos), sort_set(REF), sort_set(ALT), str(ID)]
    low_pos_GT = list()
    low_pos_REF = list()#["NA"]
    low_pos_ALT = list()#["NA"]
    for i in range(len(chrs)):
	#print i
        if chrs[i]==cur_chr and pos[i]==low_pos:
	    #print sort_set(REF), sort_set(ALT)
            low_pos_GT.append(str(GT[i])) #GT simply is a field you selected before
	    low_pos_REF.append(str(REF[i]))
	    low_pos_ALT.append(str(ALT[i]))
	else:
	    low_pos_GT.append('NA')
	    low_pos_REF.append('NA')
	    low_pos_ALT.append('NA')
	    
    info = '\t'.join([str(cur_chr),str(low_pos), sort_set(low_pos_REF), sort_set(low_pos_ALT), str(ID)])
    
    out_GT=info + '\t' + '\t'.join(low_pos_GT)
    #print sort_set(REF), lines
    outfile.write(out_GT+'\n')
    
    return None

def should_we_go(lines):
    #False when all lines are emtpy
    for i in lines:
        if i != "":
            return True
    return False    
    
def sort_set(LIST):
    while ('NA' in LIST):
	LIST.remove('NA')
    
    while (len(LIST) > 1 and '.' in LIST):
	LIST.remove('.')
    
    counter=collections.Counter(LIST)
    sort_set_list = [pair[0] for pair in sorted(counter.items(), key=lambda item: item[1], reverse=True)]
    return(','.join(sort_set_list))

if __name__=='__main__':
    import sys
    import subprocess
    import os
    import argparse
    from Bio import bgzf
    import re
    import ntpath
    import struct
    import hashlib
    import collections
    import itertools
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--input","-i"  , dest='i',required=True, action='append', help="input files")
    parser.add_argument("--headerlines","-hd"  , dest='h',required=True, type=int, help="header line number")
    parser.add_argument("--individual","-id"  , dest='id',required=True, type=str, help="Id of the individual")
    parser.add_argument("--output","-o"  , dest='o',required=True, type=str, help="output file: Matrix per sample with all tissue information")
    args = parser.parse_args()
    
    fileids=list()
    chrs=list()
    GT=list()
    REF=list()
    ALT=list()
    pos=list()
    lines=list()
    ID = args.id
    
    outfile= open(args.o,'w')
    header=['CHR','POS', 'REF', 'ALT' , 'Id']
    for id_ in range(len(args.i)):
	
        tmp=open(args.i[id_],'r')
        fileids.append(tmp)
        for i in range(args.h):
            useless = fileids[id_].readline() #skip the header
        useless = fileids[id_].readline()
        lines.append(useless)
        curline=useless.strip()
        
        GT.append(0)
	REF.append(0)
	ALT.append(0)
	
        chrs.append('x')
        pos.append(0)
	
	header.append(args.i[id_].split('/')[-1].split('.')[1])
        updatefields(GT,REF, ALT, chrs,pos,id_,curline)
	
    cur_chr=1
    cur_pos=0
    
    outfile.write('\t'.join(header)+'\n')
        
    go = should_we_go(lines) # define here an expression to describe that all lines are empty   
    
    while go:
        for i in range(len(args.i)):
	    
            while read_one_line(chrs[i],pos[i],cur_chr,cur_pos):
            #if read_one_line(chrs[i],pos[i],cur_chr,cur_pos):#check if the chr is the same and if cur_pos is greater or equal than the file
                lines[i]=fileids[i].readline().strip()
                
		if lines[i]=="":
		    chrs[i] = ''
		    pos[i] = -1
                    break
		
                GT, REF, ALT, chrs,pos = updatefields(GT, REF, ALT,chrs,pos,i,lines[i]) # updates the fields for the current line

        go = should_we_go(lines)
        #now check char
        if check_chr(cur_chr,chrs):
            
            if go:
                cur_pos = calc_low(chrs,pos,cur_chr,cur_pos) #calculate the lowest current position
                write_me(chrs,pos, REF, ALT, ID, GT, cur_chr, outfile, cur_pos)
        else:
            #print lines
            if should_we_go(lines):
              cur_chr = update_chr(chrs) #routine to update the current char, should check that's the same for all
              cur_pos = 0
            
        pass        
        
    
    
    for i in fileids:
        i.close()
    outfile.close()
