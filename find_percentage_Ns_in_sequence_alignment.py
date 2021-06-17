#!/usr/bin/env python

'''Script for finding the frequency of each nucleotide (including N and gap) in a sequence alignment'''
import sys, argparse, progressbar
from Bio import SeqIO
from time import sleep

def Usage():
	print("find_percentage_Ns_in_sequence_alignment.py -o outputfile.csv alignment.fasta")
	sys.exit(2)

def main(argv):
	parser=argparse.ArgumentParser(description="Strip sites from multi-FASTA alignments")
	parser.add_argument("alignment")
	parser.add_argument("-o", "--output_file", help="Path to output file (Default: ./output.fasta)",default="outputfile.fasta")
	
	args = parser.parse_args()
	
	print("You have supplied the following arguments:")
	print("MSA input file:\t\t" + args.alignment)
	print ("Output statistics file:\t\t" + args.output_file)
					
	# FILTER ALIGNMENT:
	with open(args.alignment, "rU") as f:
		with open(args.output_file, "w") as g:
			g.write("pos,A,C,G,T,N,-\n")
			handle = SeqIO.parse(f, "fasta")
			backup = {l.id : l.seq._data for l in handle}
			dic = backup.values()
			ids = backup.keys()

			length = len(list(dic)[0])
			num_seqs = len(ids)
			print("Alignment length: %s" % length)
			print("Number of sequences: %s" % num_seqs)
			
			print("Finished reading MSA file. Iterating over all positions and all isolates...")
			
			# Set up progress bar:
			bar = progressbar.ProgressBar(maxval=length, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
			bar.start()

			char_per_pos = {}

			for nuc in range(length): # For each nucleotide in the alignment

				col = "".join([x[nuc] for x in dic])
				
				# Make upper-case:
				col = col.replace("a","A").replace("c","C").replace("g","G").replace("t","T")
					
				nucleotide_counts = {}
				nucleotide_freqs = []
				for char in ["A","C","G","T","N","-"]:
					charcount = col.count(char)
					nucleotide_counts[char] = charcount
					freq = charcount/num_seqs
					nucleotide_freqs.append("%.3f" % freq)
				
				char_per_pos[nuc] = nucleotide_counts
				g.write(",".join([str(nuc+1)] + nucleotide_freqs) + "\n")

				bar.update(nuc+1)
			bar.finish()

if __name__ == "__main__":
	main(sys.argv[1:])
