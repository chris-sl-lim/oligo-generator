import oligo_generator.utilities.helper_functions as ogu
import itertools


class oligo_generator:

    def __init__(self, base_nt_seq='', num_changes=1):

        # Assign properties
        self.base_nt_seq = base_nt_seq
        self.num_changes = num_changes

        # Create properties
        self.generated_aa_seq = []
        self.generated_nt_seq = []
        self.generated_aa_changes = []
        self.restriction_sites = []
        self.restricted_aa_sequences = []

    def generate_aa_sequences(self, s_io = False):

        # Loop through how many changes we need
        for change in range(1, self.num_changes+1):

            # Print
            print('Generating sequences for ', change, ' number of changes.')

            # Get all combinations
            num_constant = len(self.base_aa_seq) - change
            change_list = ([True] * change) + ([False] * num_constant)
            change_combos = list(set(itertools.permutations(change_list)))

            # Create a list to store indices to remove
            idxToRemove = []

            # Combine with the change_aa_vector which says which amino acid we
            # need to change
            for idx, combo in enumerate(change_combos):
                # Create the change mask that says what elements of the 
                # change_vec needs to change
                change_mask = [a and not b for a, b in 
                               zip(list(combo), self.change_aa_vector)]
                # AND with the change_mask to apply changes
                combo = [not a and b for a, b in zip(change_mask, list(combo))]
                # Sum the amount of True values in the list
                num_true = sum(combo)
                # If the number of true values is less than the change number,
                # mark for deletion
                if num_true < change:
                    idxToRemove.append(idx)

            # Remove the ones that we need to remove
            for idx in reversed(idxToRemove):
                del change_combos[idx]

            # Now generate the amino acid sequences        
            generated_aa_seq, generated_aa_changes = \
                ogu.generate_aa_sequences(self.base_aa_seq, change_combos, self.restricted_aa_sequences)

            # Append to property
            self.generated_aa_seq = self.generated_aa_seq + generated_aa_seq
            self.generated_aa_changes = self.generated_aa_changes + \
                generated_aa_changes
            
            # Broadcast progress on SocketIO if provided.
            if s_io != False:
                progress = (change / self.num_changes) / 2 * 100
                s_io.emit('update_progress', {'progress': progress, 'current_state': change, 'total': self.num_changes})

        return

    def generate_nt_sequences(self, s_io = False):

        # Generate the nucleotide sequences
        self.generated_nt_seq, self.generated_nt_seq_change_attempts = \
            ogu.generate_nt_sequences(
                self.generated_aa_seq, self.generated_aa_changes,
                self.base_nt_seq, self.change_nt_vector,
                self.fullyfree_vector, self.restriction_sites,
                s_io = s_io
                )

        return
    

    def set_fixed_aa(self, value):

        # Loop through the sequence
        for idx, aa in enumerate(self.base_aa_seq):
            if aa == value:
                self.set_aa_pos(idx, False)


    def set_aa_pos(self, idx, value):

        # Create copy, change value, assign back
        aa_change_vec = self.change_aa_vector[:]
        aa_change_vec[idx] = value
        self.change_aa_vector = aa_change_vec[:]

        return


    def set_nt_pos(self, idx, value):

        # Create copy, change value, assign back
        nt_change_vec = self.change_nt_vector[:]
        nt_change_vec[idx] = value
        self.change_nt_vector = nt_change_vec[:]

        return

    def export_aa_sequences(self, fn):

        # Only export if the list is populated
        if len(self.generated_aa_seq) > 0:

            # Open file
            with open(fn, 'w', newline='') as file:

                # Write to csv
                # writer = csv.writer(file)
                # writer.writerow(self.generated_aa_seq)
                file.write('\n'.join(self.generated_aa_seq))

        return

    def export_nt_sequences(self, fn):

        # Only export if the list is populated
        if len(self.generated_nt_seq) > 0:

            # Open file
            with open(fn, 'w', newline='') as file:

                # Write to csv
                # writer = csv.writer(file)
                # writer.writerow(self.generated_nt_seq)
                file.write('\n'.join(self.generated_nt_seq))

        return

    @property
    def base_nt_seq(self):
        return self._base_nt_seq

    @base_nt_seq.setter
    def base_nt_seq(self, value):

        # Set the value of the base sequence
        self._base_nt_seq = value

        # Also set the base amino acid sequence
        self._base_aa_seq = ogu.nt2aa(self._base_nt_seq)

        # Set the change indices
        self._change_nt_vector = [True] * len(self.base_nt_seq)
        self._change_aa_vector = [True] * len(self.base_aa_seq)
        self._fullyfree_vector = [True] * len(self.base_aa_seq)

    @property
    def base_aa_seq(self):
        return self._base_aa_seq

    @property
    def change_nt_vector(self):
        return self._change_nt_vector

    @change_nt_vector.setter
    def change_nt_vector(self, value):
        # Assign value
        self._change_nt_vector = value
        # Sync values to the nucleotide vector
        self._change_aa_vector, self._fullyfree_vector = \
            ogu.sync_nt_change_to_aa_change(value, self._change_aa_vector)

    @property
    def change_aa_vector(self):
        return self._change_aa_vector

    @change_aa_vector.setter
    def change_aa_vector(self, value):
        # Assign value
        self._change_aa_vector = value
        # Sync values to the nucleotide vector
        self._change_nt_vector, self._fullyfree_vector = \
            ogu.sync_aa_change_to_nt_change(value, self._change_nt_vector)

    @property
    def fullyfree_vector(self):
        return self._fullyfree_vector

    @fullyfree_vector.setter
    def fullyfree_vector(self, value):
        # Assign value
        self._fullyfree_vector = value

    @property
    def restriction_sites(self):
        return self._restriction_sites

    @restriction_sites.setter
    def restriction_sites(self, value):
        self._restriction_sites = value

    @property
    def restricted_aa_sequences(self):
        return self._restricted_aa_sequences

    @restricted_aa_sequences.setter
    def restricted_aa_sequences(self, value):
        self._restricted_aa_sequences = value        
