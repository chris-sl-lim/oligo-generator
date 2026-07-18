import oligo_generator.models.generator as og
import python_codon_tables as pct
import oligo_generator.utilities.helper_functions as helper


def test_FixNucleotidePosition():

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGCTGCATT'
    o = og.oligo_generator(base_seq)

    # Get the amino acid sequence and fix position 3
    ntSeq = o.base_nt_seq
    o.set_nt_pos(5, False)
    o.set_nt_pos(8, False)

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


def test_FixedAminoAcid():

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGCTGCATT'
    o = og.oligo_generator(base_seq)

    # Fix amino acid 'C'
    o.set_fixed_aa('C')

    # Set num changes
    o.num_changes = 2

    # Generate
    o.generate_aa_sequences()
    o.generate_nt_sequences()

    # Assert that position 2 is always c
    for value in o.generated_aa_seq:
        assert(value[2] == 'C')


def test_RestrictedAminoAcidSequence():

    # Define a base sequence        
    base_seq = 'AGAAGCTGCATT'
    o = og.oligo_generator(base_seq)

    # Set num changes
    o.num_changes = 2

    # Set restricted amino acid sequence
    o.restricted_aa_sequences = ['RS']

    # Generate
    o.generate_aa_sequences()

    # Assert that 'RS' does not appear in any of the sequences
    for seq in o.generated_aa_seq:
        assert('RS' not in seq)




def test_NucleotideSequencesMatchDisplayedAminoAcidsAfterFiltering():

    # Fixed nucleotide positions can make some amino-acid candidates impossible
    # to encode. The remaining nucleotide and amino-acid rows must stay aligned.
    base_seq = 'AAAAGACTG'
    o = og.oligo_generator(base_seq)

    o.set_nt_pos(1, False)
    o.set_nt_pos(2, False)
    o.num_changes = 1

    o.generate_aa_sequences()
    o.generate_nt_sequences()

    assert len(o.generated_nt_seq) == len(o.generated_aa_seq)
    for nt_seq, aa_seq in zip(o.generated_nt_seq, o.generated_aa_seq):
        assert helper.nt2aa(nt_seq) == aa_seq


def test_LongerNucleotideGenerationWithMultipleFixedPositions():

    base_seq = 'AGA' * 6
    o = og.oligo_generator(base_seq)

    fixed_aa_positions = [2, 5]
    fixed_nt_positions = [1, 2, 7, 16, 17]

    for position in fixed_aa_positions:
        o.set_aa_pos(position, False)
    for position in fixed_nt_positions:
        o.set_nt_pos(position, False)

    o.num_changes = 3
    o.generate_aa_sequences()
    o.generate_nt_sequences()

    assert len(o.generated_nt_seq) > 0
    assert len(o.generated_nt_seq) == len(o.generated_aa_seq)
    assert len(set(o.generated_nt_seq)) == len(o.generated_nt_seq)

    fixed_aa_indices = [position - 1 for position in fixed_aa_positions]
    fixed_nt_indices = [position - 1 for position in fixed_nt_positions]

    for nt_seq, aa_seq in zip(o.generated_nt_seq, o.generated_aa_seq):
        assert len(nt_seq) == len(base_seq)
        assert helper.nt2aa(nt_seq) == aa_seq
        assert all(nt_seq[idx] == base_seq[idx] for idx in fixed_nt_indices)
        assert all(aa_seq[idx] == o.base_aa_seq[idx] for idx in fixed_aa_indices)
        assert 1 <= sum(
            aa != base_aa for aa, base_aa in zip(aa_seq, o.base_aa_seq)
        ) <= o.num_changes
