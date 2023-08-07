import python_codon_tables as pct

def nt2aa( seq ):
    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Create an empty dictionary
    codonlist = dict()

    # Create empty amino acid seq
    aa_seq = ''

    # Flatten the codon table
    for key in ct:

        # Get table for the key and loop through
        aatable = ct[key]        
        for aa in aatable:
            codonlist[aa] = key

    # Convert the sequence
    seq_len = len(seq)
    for ii in range(0, seq_len, 3):
        seq_chunk = seq[ii:ii+3]
        aa_seq = aa_seq + codonlist[seq_chunk]

    return aa_seq

def generate_aa_sequences(seq, aa_change_vec, aa_base_num_changes):
    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Get the list of available amino acids
    available_aa = list(ct.keys())

    # Create a new list for the generated sequences
    generated_aa = list()
    generated_aa_num_changes = list()

    # Remove certain amino acids
    available_aa.remove('*')
    available_aa.remove('C')
    available_aa.remove('M')

    # Loop through the elements we can change
    for idx, tf in enumerate(aa_change_vec):
        # Need to generate a new set
        if tf == 1:
            # What amino acids can we use?
            aa_to_use = available_aa
            aa_to_use.remove(seq[idx])
            # Get strings before and after
            str_before = seq[0:idx]
            str_after  = seq[idx+1:]
            # Add to the aa_num_changes vector
            aa_num_changes = aa_base_num_changes[:]
            aa_num_changes[idx] = True
            # Loop and create
            for aa in aa_to_use:
                generated_aa.append(str_before + aa + str_after)
                generated_aa_num_changes.append(aa_num_changes)

    return generated_aa, generated_aa_num_changes






