import sys
import subprocess
import os
import argparse

import annotate
import assemble


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

RESULTS_DIR = "jupiter_output/"
ANNOTATE_DIR =   RESULTS_DIR+"annos/"

JUPITER_DIR = os.path.dirname(os.path.realpath(__file__))
PRODIGAL_PROTEINS = 'predicted_proteins.faa'
PRODIGAL_GENES = 'predicted_genes.txt'

NEW_DB = JUPITER_DIR+'/db/'

def run_fastqc(forward_reads,reverse_reads,outdir):


	resultdir = outdir+'fastqc_report/'
	os.mkdir(resultdir)
	subprocess.run(["fastqc",forward_reads,reverse_reads,"-o",resultdir])


def run_spades(forward_reads,reverse_reads,outdir):


	command = ["spades.py","-1",forward_reads,"-2",reverse_reads,"--only-assembler","-o",outidr+"spades_output"]
	subprocess.run(command)

def run_trimmomatic(forward_reads,reverse_reads,forward_paired_reads,forward_unpaired_reads,reverse_paired_reads,reverse_unpaired_reads):

	command = [
		"trimmomatic","PE",
		forward_reads,reverse_reads,
		forward_paired_reads,forward_unpaired_reads,
		reverse_paired_reads,reverse_unpaired_reads,
		"LEADING:10","TRAILING:10","SLIDINGWINDOW:5:20"
	]

	subprocess.run(command)


def assemble_genome(forward_reads,reverse_reads,outdir):

	try:
		os.mkdir(outdir)

	except FileExistsError as error:
		print("Folder {} already exists.".format(outdir))
	else:

		base_Ffile = os.path.basename(forward_reads)
		base_Rfile = os.path.basename(reverse_reads)

		forward_paired = outdir+base_Ffile.replace('.','_paired.')
		forward_unpaired = outdir+base_Ffile.replace('.','_unpaired.')
		reverse_paired = outdir+base_Rfile.replace('.','_paired.')
		reverse_unpaired = outdir+base_Rfile.replace('.','_unpaired.')


	run_trimmomatic(
		forward_reads,reverse_reads,
		forward_paired,forward_unpaired,
		reverse_paired,reverse_unpaired
	)


	run_fastqc(forward_paired,reverse_paired,outdir)
		
	run_spades(forward_paired,reverse_paired,outdir)

def run_prodigal(genome,outdir):

	command = [
		'prodigal','-i',genome,'-o',outdir+PRODIGAL_GENES,
		'-a',outdir+PRODIGAL_PROTEINS
	]

	subprocess.run(command)

def makeblastdb(dbname,dbseqs):

	command = [
		'makeblastdb','-hash_index','-in',dbseqs,'-dbtype','prot',
		'-title',dbname,'-out',NEW_DB+dbname


	]


	subprocess.run(command)


def run_blastp(query,dbname,outdir):

	command = [
		"blastp","-db",dbname,"-query",query,
		'-evalue','1e-16','-max_target_seqs','2',
		'-out',outdir+"blast_results.xml","-outfmt","5"
	]

	print("Running blastp..")
	subprocess.run(command)
	print("Done")


def seq_lookup_table(fasta_file):


	lookup_table = {}
	for record in SeqIO.parse(fasta_file,"fasta"):
		lookup_table[record.id] = record.seq
	return lookup_table

def go_through(blast_record):


	prot_functions = []
	for alignment in blast_record.alignments:
		title = alignment.title
		index = title.find("sp") if "sp" in title else 0
		prot_function = title[title.find("",index):title.find('OS')]
		prot_functions.append(prot_function)
	return prot_functions


def hits_from_blast_results(result_file):

	with open(result_file) as blast_file:
		blast_records = NCBIXML.parse(blast_file)

		hits = {}
		for blast_record in blast_records:
			query = blast_record.query
			query = query[:query.find("#")].strip(" ")
			protein_functions = go_through(blast_record)
			if protein_functions:
				hits[query] = protein_functions[0]

	return hits

def label_proteins(predicted_proteins_file,blast_result_file,outfile):
	lookup_table = seq_lookup_table(predicted_proteins_file)
	hits = hits_from_blast_results(blast_result_file)

	annotations = []
	for i, item in enumerate(hits.items()):
		fasta_id,predicted_function = item
		seq = lookup_table[fasta_id]
		new_record = SeqRecord(
			id="jupiter{}".format(i),
			description="jupiter{} {}".format(i,predicted_function),
			seq=seq
	)
	annotations.append(new_record)
	SeqIO.write(annotations,outfile,"fasta")


def annotate_proteins(genome,outdir,dbname,dbseqs):

	if not os.path.isdir(NEW_DB):
		os.mkdir(NEW_DB)

	os.mkdir(outdir)
	run_prodigal(genome,outdir)
	makeblastdb(dbname,dbseqs)
	run_blastp(outdir+PRODIGAL_PROTEINS,NEW_DB+dbname,outdir)
	label_proteins(outdir+PRODIGAL_PROTEINS,
			outdir+"blast_results.xml",
			outdir+"jupiter_annotations.faa"
			)

def main():

	parser = argparse.ArgumentParser(description="Assembles and annotates genome given a set of proteins.")
	parser.add_argument("F",
				type = str,
				metavar="<forward reads>",
				help="Forward reads in FastQ format."
				)
	parser.add_argument("R",
				type=str,
				metavar="<reverse reads>",
				help="Reverse reads in FastQ format."
				)
	parser.add_argument("dbname",
				type=str,
				metavar="<database name>",
				help="Name for the blast database.")
	
	parser.add_argument("dbseqs",
				type=str,
				metavar="<database sequences>",
				help="Set of reference proteins.")


	args = parser.parse_args()

	assemble.assemble_genome(
		args.F,
		args.R,
		RESULTS_DIR
	)



	
	annotate.annotate_proteins(
		RESULTS_DIR+'spades_output/contigs.fasta',
		ANNOTATE_DIR,
		args.dbname,
		args.dbseqs
	)

if __name__ == "__main__":
	main()
