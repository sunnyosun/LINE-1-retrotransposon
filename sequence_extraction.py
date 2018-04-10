"""
This script will extract sequence from genome.fa file based on start and end positions
by Xiaoji Sun
"""

################################################################################
# Modules

# import regular expression module
import re

# import sys module
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import pandas as pd


################################################################################
# Functions
# load in hg19 or hg38 and extract only 24 chromosomes
def load_genome (genome_filename):
	# read in the genome.fa file using SeqIO
	genome=list(SeqIO.parse(genome_filename, 'fasta'))
	# only take chr1~22 + chrX + chrY
	chrs=[]
	for i in range(22):
		chrs.append('chr'+str(i+1))
	chrs.append('chrX')
	chrs.append('chrY')
	ids=[]
	for i in range(len(genome)):
		ids.append(genome[i].id)
	# put the 24 chrs into chrs_clean
	chrs_clean=[]
	for i in range(len(chrs)):
		chr=chrs[i]
		id_index=ids.index(chr)
		chrs_clean.append(genome[id_index])
	return chrs_clean

# extract sequence based on positions
def seq_extract (chrs_clean, query_filename):
	with open(query_filename) as fp:
		out_line=[]
		for i, line in enumerate(fp):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])
			end=int(line.split('\t')[2])
			index=[m for m in range(24) if chrs_clean[m].id==chr][0]
			sequence=chrs_clean[index].seq[start:end]
			strand=line.split('\t')[3]
			if strand=='-':
				sequence=sequence.reverse_complement()
			out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+sequence)
	fp.close()
	f=open(query_filename+'_seqs.fa','w+')
	f.write('\n'.join(str(a) for a in out_line))
	f.close()


# extract sequence based on positions and >= length
def seq_extract_length (chrs_clean, query_filename, length):
	with open(query_filename) as fp:
		out_line=[]
		for i, line in enumerate(fp):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])
			end=int(line.split('\t')[2])
			index=[m for m in range(24) if chrs_clean[m].id==chr][0]
			sequence=chrs_clean[index].seq[start:end]
			strand=line.split('\t')[3]
			if strand=='-':
				sequence=sequence.reverse_complement()
			if (end - start)>=length:
				out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+sequence)
	fp.close()
	print len(out_line)
	f=open(query_filename+'_'+str(length)+'_seqs.fa','w+')
	f.write('\n'.join(str(a) for a in out_line))
	f.close()


# extract L1 5'UTR sequence: before ATGGGG
def seq_extract_length_utr (chrs_clean, query_filename, length):
	with open(query_filename) as fp:
		out_line=[]
		for i, line in enumerate(fp):
			chr=line.split('\t')[0]
			start=int(line.split('\t')[1])
			end=int(line.split('\t')[2])
			index=[m for m in range(24) if chrs_clean[m].id==chr][0]
			sequence=chrs_clean[index].seq[start:end]
			strand=str(line.strip().split('\t')[3])
			print strand
			if strand == '-':
				sequence=sequence.reverse_complement()
			if (end - start)>=length:
				sequence=sequence.upper()
				pos=sequence.find('ATGGGG')
				sequence=sequence[0:pos]
				out_line.append('>'+str(chr)+':'+str(start)+'-'+str(end)+'\n'+str(sequence))
	fp.close()
	print len(out_line)
	f=open(query_filename+'_'+str(length)+'_5UTR'+'_seqs.fa','w+')
	f.write('\n'.join(str(a) for a in out_line))
	f.close()

# this function splits a fasta file into multiple files, each containing one sequence
def split_fasta (fasta_filename):
	fasta=list(SeqIO.parse(fasta_filename, 'fasta'))
	for i in range(len(fasta)):
		sequence=fasta[i].seq
		id=str('>'+fasta[i].id)
		f=open(fasta[i].id+'.fa','w+')
		f.write(id+'\n'+str(sequence))
		f.close()



##########################################################################


#genome_filename='hg19.genome.fa'
genome_filename='hg38.fa'
chrs_clean=load_genome(genome_filename)

query_filename='ucsc.rmsk.hg38.L1HS.bed'
length=5000
seq_extract_length_utr(chrs_clean, query_filename, length)








