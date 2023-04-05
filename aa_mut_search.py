#!/usr/bin/python3

"""
Program to find all the single nucleotide changes possible between two
amino acid changes given a starting codon
"""

import sys
import argparse
import re

# Set codon dictionary
codon_dict = {'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L',
'CTA':'L', 'CTG':'L', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M', 'GTT':'V',
'GTC':'V', 'GTA':'V', 'GTG':'V', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'ACT':'T', 'ACC':'T', 'ACA':'T',
'ACG':'T', 'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 'TAT':'Y', 'TAC':'Y',
'TAA':'X', 'TAG':'X', 'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q', 'AAT':'N',
'AAC':'N', 'AAA':'K', 'AAG':'K', 'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
'TGT':'C', 'TGC':'C', 'TGA':'X', 'TGG':'W', 'CGT':'R', 'CGC':'R', 'CGA':'R',
'CGG':'R', 'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R', 'GGT':'G', 'GGC':'C',
'GGA':'G', 'GGG':'G'}

def run(args):
	codon = args.codon.upper()
	variant = args.variant.upper()
	wt = codon_dict[codon]
	index = 0
	while index <= 2:
		for mutation in ['A', 'G', 'C', 'T']:
			if index == 0:
				mutated_codon = re.sub("(\w)(\w)(\w)", f"{mutation}"+r"\2\3", codon)
				mutated_aa = codon_dict[mutated_codon]
				if(mutated_aa == variant):
					print(f"{codon}>{mutated_codon} = {wt}>{mutated_aa}")
			if index == 1:
				mutated_codon = re.sub("(\w)(\w)(\w)", r"\1"+f"{mutation}"+r"\3", codon)
				mutated_aa = codon_dict[mutated_codon]
				if(mutated_aa == variant):
					print(f"{codon}>{mutated_codon} = {wt}>{mutated_aa}")
			if index == 2:
				mutated_codon = re.sub("(\w)(\w)(\w)", r"\1\2"+f"{mutation}",  codon)
				mutated_aa = codon_dict[mutated_codon]
				if(mutated_aa == variant):
					print(f"{codon}>{mutated_codon} = {wt}>{mutated_aa}")
		index += 1

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-c',dest='codon', help='The inpurt codon that will be mutated', 
	required=True)
	parser.add_argument('-v',dest='variant', help='The mutated amino acid',
	required=True)
	args = parser.parse_args()
	run(args)

