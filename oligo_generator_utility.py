import python_codon_tables as pct
import fnmatch

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

def generate_aa_sequences(seq, nt_seq, aa_change_vec, nt_change_vec, aa_base_num_changes):
    # seq: base amino acid sequence
    # nt_seq: base nucleotide sequence
    # aa_change_vec: boolean the same size as the sequence which tells us what amino acids we can change
    # nt_change_vec: boolean the same size as the nucleotide sequence which tells us what nucleotides we can change
    # aa_base_num_changes: boolean the same size as the sequence which tells us what amino acids have already changed

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
            # Get corresponding nt sequence
            nt_change_chunk = nt_change_vec[(idx*3):(idx*3)+3]

            # What amino acids can we use? If all the nt_chunk is True, then we are fully free.
            # If not, then we need to get a sub of available amino acids
            if nt_change_chunk == [True]*3:
                aa_to_use = available_aa
                
            else:
                nt_seq_chunk = seq[(idx*3):(idx*3)+3]
                aa_to_use = get_compatible_aa(nt_seq_chunk, nt_change_chunk)

            # Remove the amino acid that is already there
            aa_to_use.remove(seq[idx])

            # Get strings before and after
            str_before = seq[0:idx]
            str_after  = seq[idx+1:]

            # Add to the aa_num_changes vector by creating a copy of the original
            aa_num_changes = aa_base_num_changes[:]
            aa_num_changes[idx] = True

            # Loop and create
            for aa in aa_to_use:
                generated_aa.append(str_before + aa + str_after)
                generated_aa_num_changes.append(aa_num_changes)

    return generated_aa, generated_aa_num_changes

def sync_nt_change_to_aa_change(nt_change_vec, aa_change_vec):

        # Create fully free
    fully_free_vec = [False] * len(aa_change_vec)

    # Loop through amino acid change vector
    for idx, aa_change in enumerate(aa_change_vec):

        # Get the corresponding nucleotide changes
        nt_change = nt_change_vec[(idx*3):(idx*3)+3]

        if nt_change == [False] * 3:

            # The corresponding amino acid should be set to Talse
            aa_change_vec[idx] = False
            fully_free_vec[idx] = False

        elif nt_change == [True] * 3:

            # The corresponding amino acid should be set to True
            aa_change_vec[idx] = True
            fully_free_vec[idx] = True

        else:

            # Change vector should be True, fully free as False
            aa_change_vec[idx] = True
            fully_free_vec[idx] = False

    return aa_change_vec, fully_free_vec

def sync_aa_change_to_nt_change(aa_change_vec, nt_change_vec):

    # Create fully free
    fully_free_vec = [False] * len(aa_change_vec)

    # Loop through amino acid change vector
    for idx, aa_change in enumerate(aa_change_vec):

        # Get the corresponding nucleotide changes
        nt_change = nt_change_vec[(idx*3):(idx*3)+3]

        if aa_change == False:
        
            # All should be set to false (if we cannot change the amino acid, we cannot change the nucleotide)
            nt_change_vec[(idx*3):(idx*3)+3] = [False] * 3

        elif aa_change == True:

            # If they are all set to false, then set them all to true.
            # If they are not all false, then they're either set correctly (e.g. 111, 101, 011) as at least one nucleotide can change
            if nt_change == [False] * 3:                
                nt_change_vec[(idx*3):(idx*3)+3] = [True]*3
                fully_free_vec[idx] = True

            elif nt_change == [True] * 3:
                fully_free_vec[idx] = True

    return nt_change_vec, fully_free_vec


def get_compatible_aa(nt_seq, nt_change):

    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Create an empty dictionary
    codonlist = list()
    codonseq = list()

    # Flatten the codon table
    for key in ct:

        # Get table for the key and loop through
        aatable = ct[key]        
        for aa in aatable:
            codonlist.append(key)
            codonseq.append(aa)

    # Convert the sequence to wildcards
    nt_search = ''
    for idx, change in enumerate(nt_change):
        if change == False:
            nt_search = nt_search + '*'
        else:
            nt_search = nt_search + nt_seq[idx]

    # Filter the codon table
    matchidx = [i for i, j in enumerate(codonseq) if j in fnmatch.filter(codonseq, nt_search)]

    # Return compatible codons
    available_aa = [codonlist[i] for i in matchidx]

    # Now get unique values
    return list(set(available_aa))




