import oligo_generator as og
import oligo_generator_utility as ogu

# define base sequence
base_seq = 'AGAAGCTGCATTGATACCATCCCCAAAAGCAGATGCACCGCCTTCCAGTGCGCGCACAGCATGAAATACCGGCTGAGCTTCTGTAGAAAGACCTGCGGCACCTGT'

# create OligoGenerator object
o = og.oligo_generator( base_seq )

# print
print(o.base_aa_seq)
print(o.base_nt_seq)

# Generate a new change vector
o.set_aa_pos(4, True)
o.set_aa_pos(5, True)
o.set_nt_pos(2, True)
o.set_nt_pos(5, True)
o.num_changes = 2

# create new amino acid sequences
o.generate_aa_sequences()

# Print generated sequences
print(o.generated_aa_seq)