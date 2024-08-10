#!/usr/bin/env python3

from Bio import SeqIO

def filter_function(): 

    parser = SeqIO.parse("$fastq_file", "fastq")

    with open("$fastq_name"+"_filtered.fastq","w") as out_handle: 
    	for record_fq in parser: 
    		if len(record_fq.seq) > 300 and len(record_fq.seq) < 800: 
        		SeqIO.write(record_fq, out_handle,"fastq") 
    return(out_handle)

filter_function() 