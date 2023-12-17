import python_codon_tables as pct
import fnmatch

from itertools import combinations


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

    # Create a new list for the generated sequences
    generated_aa = list()
    generated_aa_changes = list()

    # For each combo
    for combo in change_combos:

        # Create a list for the static chunks (e.g. bits of the sequence that
        # do not change)
        static_chunks = list()
        chunk = ''

        # Create a list for an amino acid library (e.g. what amino acids can
        # be used at what position)
        aa_lib = list()

        # Build up a list of the static chunks
        for idx, tf in enumerate(combo):

            if tf is False:
                chunk = chunk + aa_seq[idx]
            elif tf is True:
                # Append the chunk to the library and reset
                static_chunks.append(chunk)
                chunk = ''
                # Take stock of what position we are at and remove the amino 
                # acid that is already there
                aa_to_use = available_aa
                if aa_seq[idx] in aa_to_use:
                    aa_to_use.remove(aa_seq[idx])
                aa_lib.append(aa_to_use)

        # Add the last chunk after the loop has finished
        static_chunks.append(chunk)

        # Create a list of roots
        prev_seq_roots = list()

        # Now create the sequences
        for idx, chunk in enumerate(static_chunks):

            # Append the latest chunk to all the generated sequence roots
            # so far
            if idx == 0:
                prev_seq_roots.append(chunk)
            else:
                for idx1, seq_root in enumerate(prev_seq_roots):
                    prev_seq_roots[idx1] = seq_root + chunk

            # If we're not at the last segment, then append the amino acids
            # as well
            if idx < (len(static_chunks)-1):

                # Get the aa_lib
                aa_to_use = aa_lib[idx]

                # Create a new_seq_roots list
                new_seq_roots = list()

                # Append
                for seq_root in prev_seq_roots:

                    # Generate a new sequence root for each amino acid
                    for aa in aa_to_use:
                        new_seq_roots.append(seq_root + aa)

                # Replace prev_seq_roots
                prev_seq_roots = new_seq_roots

            else:

                # Now we're at the end - these are our sequences
                generated_aa = generated_aa + prev_seq_roots
                for ii in range(0, len(prev_seq_roots)):
                    generated_aa_changes.append(list(combo))

        # Now check for restricted amino acid sequences
        idx_to_remove = []
        for restricted_seq in restricted_aa_seq:
            for idx, seq in enumerate(generated_aa):
                if restricted_seq in seq:
                    idx_to_remove.append(idx)

        # Recreate list
        generated_aa = [j for i, j in enumerate(generated_aa) if i not in idx_to_remove]

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
                          s_io = False):

    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Create list for generated nucleotide sequences
    generated_nt = list()
    generated_change_attempts = list()

    # We have to loop through each sequence
    for idx, aa_seq in enumerate(aa_sequences):

        # Create a vector to track how many options we have tried in each
        # position
        aa_num_change_attempts = [0] * len(aa_seq)

        # Set valid_seq flag
        valid_seq = False

        # Enter while loop
        while valid_seq is False:

            # Create a blank nucleotide sequence
            nt_seq = ''

            # Get the vector that describes what amino acids have changed
            aa_change = aa_num_changes[idx]

            # Now we have to through each aa in the sequence
            for idx_pos, aa in enumerate(aa_seq):

                if aa_change[idx_pos] is False:

                    # If we have set the amino acid to stay the same - copy
                    # the chunk of the base nt seq across
                    nt_seq = nt_seq + base_nt_seq[(idx_pos*3):(idx_pos*3)+3]

                elif aa_change[idx_pos] is True:

                    # Get the bit of the codon table pertaining to that amino
                    # acid and sort by probability
                    ct_slice = sorted(ct[aa].items(), key=lambda item:
                                      item[1], reverse=True)

                    # If we have changed the amino acid -- are we fully free?
                    if fullyfree_vec[idx_pos] is True:

                        # Check whether we have gone through all the 
                        # possibilities
                        if aa_num_change_attempts[idx_pos] >= len(ct_slice):
                            nt_seq = ''
                            break

                        # Insert the most prevalant nt translation of that aa
                        # into the sequence
                        nt_chunk = ct_slice[aa_num_change_attempts[idx_pos]][0]
                        nt_seq = nt_seq + nt_chunk

                        # Increment aa_num_changes vector
                        aa_num_change_attempts[idx_pos] = \
                            aa_num_change_attempts[idx_pos] + 1

                    else:

                        # If we are not fully free, get the base nt sequence
                        base_nt_chunk = base_nt_seq[(idx_pos*3):(idx_pos*3)+3]

                        # Get what we are able to change
                        nt_change_chunk = \
                            nt_change_vec[(idx_pos*3):(idx_pos*3)+3]

                        # Get a side by side list
                        nt_list = []
                        for item in ct_slice:
                            nt_list.append(item[0])

                        # Filter the codon list
                        matchidx = filter_codon_list(base_nt_chunk,
                                                     nt_change_chunk, nt_list)

                        # Get the subset of the list
                        ct_subset = []
                        for ii in matchidx:
                            ct_subset.append(nt_list[ii])

                        # Check whether we have gone through all the
                        # possibilities
                        if aa_num_change_attempts[idx_pos] >= len(ct_subset):
                            nt_seq = ''
                            break
                        # Insert the most prevalant nt translation of that aa
                        # into the sequence
                        nt_chunk = ct_subset[aa_num_change_attempts[idx_pos]]
                        nt_seq = nt_seq + nt_chunk
                        # Increment aa_num_changes vector
                        aa_num_change_attempts[idx_pos] = \
                            aa_num_change_attempts[idx_pos] + 1

            # Set valid_seq flag
            valid_seq = True

            # If we don't have any sequence, then bail out, can't find a
            # solution so we don't need to look for restriction sites
            if nt_seq != '':

                # Now we need to check for restricton sites
                for res_site in restriction_sites:

                    # if restriction site exists in the sequence, set flag
                    if res_site in nt_seq:
                        valid_seq = False

            # Print message
            print('Generated nucleotide sequence ', idx+1, ' of ',
                  len(aa_sequences))
            
            # Broadcast progress on SocketIO if provided.
            if (s_io != False and (((idx+1) % 10 == 0)) or idx+1 == len(aa_sequences)):
                progress = 50 + (((idx+1) / len(aa_sequences)) / 2 * 100)
                s_io.emit('update_progress', {'progress': progress, 'current_state': idx+1, 'total': len(aa_sequences)})

        # Exited the while loop so store the sequence
        generated_nt.append(nt_seq)
        generated_change_attempts.append(aa_num_change_attempts)

    # Remove empties from list
    nonemptyidx = [generated_nt.index(x) for x in generated_nt if x != '']
    generated_nt = [generated_nt[i] for i in nonemptyidx]
    generated_change_attempts = [generated_change_attempts[i] for i in nonemptyidx]


    # Gone through all the amino acid sequences
    return generated_nt, generated_change_attempts
