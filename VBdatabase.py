import sys
import argparse


from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():

	parser = argparse.ArgumentParser(descrip='Get the virulence genes at http://www.mgc.ac.cn/VFs/Down/CP_VFs.ffn.gz)

	parser.add_argument('--infile',

						required = True,
	parser.add_argument('--genus',
						required = False,

	return parser.parse_args()

def main():
	
	args = parse_args()

	Vb = {}


	for record in SeqIO.parse(open(args.infile, "r"), "fasta"):
		name = record.description
		genus = name.split("[")[-1].split()[0]
		if (not args.genus) or (genus == args.genus):

			if genus in Vb:
				Vb[genus].append(record)

			else:
				 Vb[genus] = [record]


	for genus in Vb:
		records = Vb[genus]
		SeqIO.write(records, (genus + ".fsa"), "fasta")



	if __name__ == '_main_':

		main()
