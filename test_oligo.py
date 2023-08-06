import oligo_generator as og

# define base sequence
base_seq = 'AGAAGCTGCATTGATACCATCCCCAAAAGCAGATGCACCGCCTTCCAGTGCGCGCACAGCATGAAATACCGGCTGAGCTTCTGTAGAAAGACCTGCGGCACCTGT'

# create OligoGenerator object
o = og.oligo_generator( base_seq)

# print
print(o.base_aa_seq)
print(o.base_nt_seq)

# Generate a new change vector
o.change_aa_vector[4] = 1

# create new amino acid sequences
o.generate_aa_sequences()

# Print generated sequences
print(o.generated_aa_seq)