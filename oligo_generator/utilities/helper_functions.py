import python_codon_tables as pct
import fnmatch

from itertools import combinations, product


def nt2aa(seq):
    """
    Converts nucleotide sequence to an amino acid sequence
    """

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


def aa2nt(seq):
    """
    Converts amino acid sequences to nucleotide sequences using highest
    likelihood nucleotide sequence
    """

    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Create blank nt_seq
    nt_seq = ''

    # Return the nt seq
    for aa in seq:

        # Get the dictionary for that aa
        nt_dict = ct[aa]
        # Sort by most likely
        nt_dict = dict(sorted(nt_dict.items(),
                              key=lambda x: x[1], reverse=True))
        # Append most likely
        nt_seq = nt_seq + list(nt_dict.items())[0][0]

    return nt_seq


def all_combinations(items_list):

    result = []

    for i in range(1, len(items_list) + 1):
        for combo in combinations(items_list, i):
            result.append(combo)

    return result


def generate_aa_sequences(aa_seq, change_combos, restricted_aa_seq):
    """
    seq: base amino acid sequence
    nt_seq: base nucleotide sequence
    aa_change_vec: boolean the same size as the sequence which tells us what
        amino acids we can change
    nt_change_vec: boolean the same size as the nucleotide sequence which
        tells us what nucleotides we can change
    aa_base_num_changes: boolean the same size as the sequence which tells us
        what amino acids have already changed
    """

    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Get the list of available amino acids
    available_aa = list(ct.keys())

    # Remove the amino acids we will never want
    if '*' in available_aa:
        available_aa.remove('*')
    if 'C' in available_aa:
        available_aa.remove('C')
    if 'M' in available_aa:
        available_aa.remove('M')

    generated_aa = []
    generated_aa_changes = []
    restricted_aa_seq = [seq for seq in restricted_aa_seq if seq]

    for combo in change_combos:
        aa_choices = []

        for idx, should_change in enumerate(combo):
            if should_change is True:
                aa_choices.append([aa for aa in available_aa
                                   if aa != aa_seq[idx]])
            else:
                aa_choices.append([aa_seq[idx]])

        for candidate in product(*aa_choices):
            candidate_seq = ''.join(candidate)
            if any(restricted_seq in candidate_seq
                   for restricted_seq in restricted_aa_seq):
                continue

            generated_aa.append(candidate_seq)
            generated_aa_changes.append(list(combo))

    return generated_aa, generated_aa_changes

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

        if aa_change is False:

            # All should be set to false (if we cannot change the amino acid,
            # we cannot change the nucleotide)
            nt_change_vec[(idx*3):(idx*3)+3] = [False] * 3

        elif aa_change is True:

            # If they are all set to false, then set them all to true.
            # If they are not all false, then they're either set correctly
            # (e.g. 111, 101, 011) as at least one nucleotide can change
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

    # Filter the codon list
    matchidx = filter_codon_list(nt_seq, nt_change, codonseq)

    # Return compatible codons
    available_aa = [codonlist[i] for i in matchidx]

    # Now get unique values
    return list(set(available_aa))


def filter_codon_list(nt_seq, nt_change, nt_list):

    # Convert the sequence to wildcards
    nt_search = ''
    for idx, change in enumerate(nt_change):
        if change is True:
            nt_search = nt_search + '?'
        else:
            nt_search = nt_search + nt_seq[idx]

    # Filter the codon table
    matchidx = [i for i, j in enumerate(nt_list) if j in
                fnmatch.filter(nt_list, nt_search)]

    return matchidx


def generate_nt_sequences(aa_sequences, aa_num_changes, base_nt_seq,
                          nt_change_vec, fullyfree_vec, restriction_sites, 
                          s_io=None):

    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    generated_nt = []
    generated_change_attempts = []
    generated_aa_indices = []
    restriction_sites = [site for site in restriction_sites if site]

    # We have to loop through each sequence
    for idx, aa_seq in enumerate(aa_sequences):

        codon_options = []
        aa_change = aa_num_changes[idx]

        for idx_pos, aa in enumerate(aa_seq):
            if aa_change[idx_pos] is False:
                codon_options.append(
                    [base_nt_seq[(idx_pos*3):(idx_pos*3)+3]]
                )
                continue

            ct_slice = sorted(ct[aa].items(), key=lambda item: item[1],
                              reverse=True)

            if fullyfree_vec[idx_pos] is True:
                codon_options.append([item[0] for item in ct_slice])
                continue

            base_nt_chunk = base_nt_seq[(idx_pos*3):(idx_pos*3)+3]
            nt_change_chunk = nt_change_vec[(idx_pos*3):(idx_pos*3)+3]
            nt_list = [item[0] for item in ct_slice]
            matchidx = filter_codon_list(base_nt_chunk, nt_change_chunk,
                                         nt_list)
            codon_options.append([nt_list[ii] for ii in matchidx])

        nt_seq = ''
        aa_num_change_attempts = [0] * len(aa_seq)

        if all(codon_options):
            indexed_options = [
                list(enumerate(options)) for options in codon_options
            ]

            for candidate in product(*indexed_options):
                nt_seq = ''.join(codon for _, codon in candidate)

                if any(res_site in nt_seq for res_site in restriction_sites):
                    continue

                aa_num_change_attempts = [
                    option_idx for option_idx, _ in candidate
                ]
                break
            else:
                nt_seq = ''

        # Broadcast progress on SocketIO if provided.
        if s_io is not None:
            if ((idx+1) % 10 == 0) or (idx+1 == len(aa_sequences)):
                progress = 50 + (((idx+1) / len(aa_sequences)) / 2 * 100)
                s_io.emit('update_progress', {'progress': progress,
                                              'current_state': idx+1,
                                              'total': len(aa_sequences)})

        if nt_seq:
            generated_nt.append(nt_seq)
            generated_change_attempts.append(aa_num_change_attempts)
            generated_aa_indices.append(idx)

    # Gone through all the amino acid sequences
    return generated_nt, generated_change_attempts, generated_aa_indices
