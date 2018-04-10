"""
This script compares the sequence similarity between L1 consensus sequences
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
import regex

################################################################################
# Functions
def nbp_scanning (query_filename, subject_filename, nbp):
	n=nbp
	q=list(SeqIO.parse(query_filename, 'fasta'))[0]
	q=q.upper()
	s=list(SeqIO.parse(subject_filename, 'fasta'))[0]
	s=s.upper()
	out=[]
	for i in range(len(q)-nbp+5):
		tmp=str(q.seq[i:i+nbp])
		pos=i+5
		f0=regex.findall(tmp, str(s.seq))
		f1=regex.findall('('+tmp+')'+'{e<=1}', str(s.seq))
		if (len(f0)>=1):
			score=1
		elif (len(f0)==0 and len(f1)>=1):
			score=0.5
		else:
			score=0
		out.append(str(pos)+'\t'+str(score))
	f=open(query_filename+'_'+subject_filename+'.txt','w+')
	f.write('\n'.join(str(a) for a in out))
	f.close()


#########################################################################
nbp_scanning('L1HS.fa','LAPA1.fa',10)






