import oligo_generator.models.generator as og
import python_codon_tables as pct


def test_FixNucleotidePosition():

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGCTGCATT'
    o = og.oligo_generator(base_seq)

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


def test_ThreeChangeFourCodonNucleotideGeneration():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to
    # get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGAAGAAGA'
    o = og.oligo_generator(base_seq)

    # Set number of changes
    o.num_changes = 3
    o.generate_aa_sequences()

    # Generate nucleotide sequences
    o.generate_nt_sequences()

    # Assert that the number of nucleotide sequences is equal to amino acid sequences
    assert len(o.generated_nt_seq) == len(o.generated_aa_seq)

    # Assert that all nucleotide sequences are unique
    assert len(set(o.generated_nt_seq)) == len(o.generated_nt_seq)