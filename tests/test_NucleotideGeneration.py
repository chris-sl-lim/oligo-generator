import oligo_generator.oligo_generator as og
import csv

from math import comb

def test_FixNucleotidePosition():
    
    # Define base seq and create oligo generator object
    base_seq = 'AGAAGCTGCATT'
    o = og.oligo_generator( base_seq )

    # Get the amino acid sequence and fix position 3
    ntSeq = o.base_nt_seq
    o.set_nt_pos(4, False)
    o.set_nt_pos(7, False)

    # Generate two changes
    o.num_changes = 2
    o.generate_aa_sequences()

    # Generate nucleotide sequences
    o.generate_nt_sequences()

    # Assert all the AA sequences are unique
    unique_seq = set(o.generated_nt_seq)
    assert len(unique_seq) == len(o.generated_nt_seq)

    # Assert that the correct position is fixed
    for seq in o.generated_nt_seq:
        if seq != '':
            assert seq[4] == ntSeq[4]
            assert seq[7] == ntSeq[7]