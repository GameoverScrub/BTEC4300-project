import re
import sys
import argparse

from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def parse_args():
	parser = argparse.ArgumentParser(description='Virulent gene extractor using the VFBD database')

	parser.add_argument('--file',
			required = True,
			help = 'Extract pure virulent genes')
	parser.add_argument('--genus',
			required = False,
			help = 'Extract the genes using the genus')
	return parser.parse_args()  

def main():
	args = parse_args()
	
	vb = {} 
	
	for record in SeqIO.parse(open(args.file, "r"), "fasta"):
		viral_name = record.description
		genus = viral_name.split("[")[-1].split()[0]
		if (not args.genus) or (genus == args.genus):
			if genus in vb:
				vb[genus].append(record)
			else:
				vb[genus] = [record]

	
	for genus in vb:
		records = vb[genus] 
		SeqIO.write(records, (genus + ".faa"), "fasta")

if __name__ == '__main__':
	main()
