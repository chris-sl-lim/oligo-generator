import oligo_generator as og
import python_codon_tables as pct

from math import comb

def test_FixAminoAcidPosition():
    
    # Define base seq and create oligo generator object
    base_seq = 'AGAAGCTGCATT'
    o = og.oligo_generator( base_seq )

    # Get the amino acid sequence and fix position 3
    aaSeq = o.base_aa_seq
    pos = 2
    o.set_aa_pos(pos, False)

    # Generate two changes
    o.num_changes = 2
    o.generate_aa_sequences()

    # Assert all the AA sequences are unique
    unique_seq = set(o.generated_aa_seq)
    assert len(unique_seq) == len(o.generated_aa_seq)

    # Assert that the correct position is fixed
    for seq in o.generated_aa_seq:
        assert seq[pos] == aaSeq[pos]

def test_OneChangeOneCodon():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGA'
    o = og.oligo_generator( base_seq )

    # Set number of changes
    o.num_changes = 1
    o.generate_aa_sequences()

    # Get number of changes
    assert len(o.generated_aa_seq) == num_ct

def test_OneChangeTwoCodon():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGA'
    o = og.oligo_generator( base_seq )

    # Set number of changes
    o.num_changes = 1
    o.generate_aa_sequences()

    # Get number of changes
    assert len(o.generated_aa_seq) == num_ct * 2   

def test_OneChangeThreeCodon():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGAAGA'
    o = og.oligo_generator( base_seq )

    # Set number of changes
    o.num_changes = 1
    o.generate_aa_sequences()

    # Get number of changes
    assert len(o.generated_aa_seq) == num_ct * 3   

def test_TwoChangeTwoCodon():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGA'
    o = og.oligo_generator( base_seq )

    # Set number of changes
    o.num_changes = 2
    o.generate_aa_sequences()

    # Get number of changes
    assert len(o.generated_aa_seq) == (comb(2,1) * num_ct) + (comb(2,2) * num_ct**2)

def test_TwoChangeThreeCodon():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGAAGA'
    o = og.oligo_generator( base_seq )

    # Set number of changes
    o.num_changes = 2
    o.generate_aa_sequences()

    # Get number of changes
    assert len(o.generated_aa_seq) == (comb(3,1) * num_ct) + (comb(3,2) * num_ct**2)

def test_ThreeChangeFourCodon():

    # Get the codon list and how many 
    ct = pct.get_codons_table("h_sapiens_9606")
    num_ct = len(ct)

    # Subtract the codons that we remove (*, C, M and the original codon) to get the num changes we expect
    num_ct = num_ct - 4

    # Define base seq and create oligo generator object
    base_seq = 'AGAAGAAGAAGA'
    o = og.oligo_generator( base_seq )

    # Set number of changes
    o.num_changes = 3
    o.generate_aa_sequences()

    # Get number of changes
    assert len(o.generated_aa_seq) == (comb(4,1) * num_ct) + (comb(4,2) * num_ct**2) + (comb(4,3) * num_ct**3)


    

    

    
