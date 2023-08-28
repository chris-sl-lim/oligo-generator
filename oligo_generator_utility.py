import python_codon_tables as pct
import fnmatch

from copy import deepcopy

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

    # Loop through the elements we can change
    for idx, tf in enumerate(aa_change_vec):

        # Need to generate a new set
        if tf == True:
            # Get corresponding nt sequence
            nt_change_chunk = nt_change_vec[(idx*3):(idx*3)+3]

            # What amino acids can we use? If all the nt_chunk is True, then we are fully free.
            # If not, then we need to get a sub of available amino acids
            if nt_change_chunk == [True]*3:
                aa_to_use = available_aa
                
            else:
                nt_seq_chunk = nt_seq[(idx*3):(idx*3)+3]
                aa_to_use = get_compatible_aa(nt_seq_chunk, nt_change_chunk)

            # Remove the amino acid that is already there
            if seq[idx] in aa_to_use:
                aa_to_use.remove(seq[idx])

            # Remove certain amino acids
            if '*' in aa_to_use:
                aa_to_use.remove('*')
            if 'C' in aa_to_use:
                aa_to_use.remove('C')
            if 'M' in aa_to_use:
                aa_to_use.remove('M')

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
        if change == False:
            nt_search = nt_search + '?'
        else:
            nt_search = nt_search + nt_seq[idx]

    # Filter the codon table
    matchidx = [i for i, j in enumerate(nt_list) if j in fnmatch.filter(nt_list, nt_search)]

    return matchidx

def generate_nt_sequences(aa_sequences, aa_num_changes, base_nt_seq, nt_change_vec, fullyfree_vec, restriction_sites):

    # Get the codon table for 'h_sapiens_9606"
    ct = pct.get_codons_table("h_sapiens_9606")

    # Create list for generated nucleotide sequences
    generated_nt = list()
    generated_change_attempts = list()

    # We have to loop through each sequence
    for idx, aa_seq in enumerate(aa_sequences):

        # Create a vector to track how many options we have tried in each position
        aa_num_change_attempts = [0] * len(aa_seq)

        # Set valid_seq flag
        valid_seq = False

        # Enter while loop
        while valid_seq == False:

            # Create a blank nucleotide sequence
            nt_seq = ''
            
            # Get the vector that describes what amino acids have changed
            aa_change = aa_num_changes[idx]

            # Now we have to through each aa in the sequence
            for idx_pos, aa in enumerate(aa_seq):

                if aa_change[idx_pos] == False:

                    # If we have set the amino acid to stay the same - copy the chunk of the base nt seq across
                    nt_seq = nt_seq + base_nt_seq[(idx_pos*3):(idx_pos*3)+3]

                elif aa_change[idx_pos] == True:

                    # Get the bit of the codon table pertaining to that amino acid and sort by probability
                    ct_slice = sorted(ct[aa].items(), key=lambda item: item[1], reverse=True)

                    # If we have changed the amino acid -- are we fully free?
                    if fullyfree_vec[idx_pos] == True:

                        # Check whether we have gone through all the possibilities
                        if aa_num_change_attempts[idx_pos] >= len(ct_slice):
                            nt_seq = ''
                            break

                        # Insert the most prevalant nt translation of that aa into the sequence
                        nt_chunk = ct_slice[aa_num_change_attempts[idx_pos]][0]
                        nt_seq = nt_seq + nt_chunk

                        # Increment aa_num_changes vector
                        aa_num_change_attempts[idx_pos] = aa_num_change_attempts[idx_pos] + 1

                    else:

                        # If we are not fully free, get the base nt sequence
                        base_nt_chunk = base_nt_seq[(idx_pos*3):(idx_pos*3)+3]

                        # Get what we are able to change
                        nt_change_chunk = nt_change_vec[(idx_pos*3):(idx_pos*3)+3]

                        # Get a side by side list
                        nt_list = []
                        for item in ct_slice:
                            nt_list.append(item[0])

                        # Filter the codon list
                        matchidx = filter_codon_list(base_nt_chunk, nt_change_chunk, nt_list)

                        # Get the subset of the list
                        ct_subset = []
                        for ii in matchidx:
                            ct_subset.append(nt_list[ii])

                        # Check whether we have gone through all the possibilities
                        if aa_num_change_attempts[idx_pos] >= len(ct_subset):
                            nt_seq = ''
                            break
                        # Insert the most prevalant nt translation of that aa into the sequence
                        nt_chunk = ct_subset[aa_num_change_attempts[idx_pos]]
                        nt_seq = nt_seq + nt_chunk
                        # Increment aa_num_changes vector
                        aa_num_change_attempts[idx_pos] = aa_num_change_attempts[idx_pos] + 1

            # Set valid_seq flag
            valid_seq = True                    

            # If we don't have any sequence, then bail out, can't find a solution so 
            # we don't need to look for restriction sites
            if nt_seq != '':

                # Now we need to check for restricton sites
                for res_site in restriction_sites:

                    # if restriction site exists in the sequence, set flag
                    if res_site in nt_seq:
                        valid_seq = False

            # Print message
            print('Generated nucleotide sequence ', idx+1, ' of ', len(aa_sequences))

        # Exited the while loop so store the sequence
        generated_nt.append(nt_seq)
        generated_change_attempts.append(aa_num_change_attempts)

    # Gone through all the amino acid sequences
    return generated_nt, generated_change_attempts






