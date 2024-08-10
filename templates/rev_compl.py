#!/usr/bin/env python3

from Bio import SeqIO

def reverse_compl(): 
    parser = SeqIO.parse("$rev_file", "fastq")
    with open("$fastq_name"+"_rev_fov.fq","w") as out_handle: 
        for record_fq in parser: 
            rev_rec = record_fq.reverse_complement() 
            SeqIO.write(rev_rec, out_handle,"fastq") 
    return(out_handle)
reverse_compl() 
