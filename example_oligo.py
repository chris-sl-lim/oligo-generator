import oligo_generator as og
from oligo_generator import oligo_generator_utility as ogu
import csv

# define base sequence
base_seq = 'AGAAGCTGCATTGATACCATCCCCAAAAGCAGATGCACCGCCTTCCAGTGCGCGCACAGCATGAAATACCGGCTGAGCTTCTGTAGAAAGACCTGCGGCACCTGT'

# create OligoGenerator object
o = og.generator( base_seq )

# print
print(o.base_aa_seq)
print(o.base_nt_seq)

# Fix certain position in the amino acid and nucleotide chains
o.set_aa_pos(4, False)
o.set_aa_pos(5, False)
o.set_nt_pos(3, False)
o.set_nt_pos(5, False)
o.num_changes = 2

# create new amino acid sequences
o.generate_aa_sequences()

# import restriction sites
file = open('RestrictionSites.csv')
csvreader = csv.reader(file)
restriction_sites = []
for row in csvreader:
    restriction_sites.append(row[1])

# set restriction_sites in the generator
o.restriction_sites = restriction_sites

# Generate nucleotide sequences
o.generate_nt_sequences()

# Print generated sequences
print(o.generated_aa_seq)