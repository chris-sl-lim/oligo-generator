# import statements
import python_codon_tables as pct
import oligo_generator_utility as ogu
import csv

class oligo_generator:

    def __init__(self, base_nt_seq='', num_changes=1):
        
        # Assign properties
        self.base_nt_seq = base_nt_seq
        self.num_changes = num_changes

        # Create properties
        self.generated_aa_seq = []
        self.generated_nt_seq = []

    def generate_aa_sequences(self):

        # Create the first array of new amino acid sequences
        self.generated_aa_seq, self.generated_aa_num_changes = ogu.generate_aa_sequences(self.base_aa_seq, self.base_nt_seq, self.change_aa_vector, self.change_nt_vector, [False]*len(self.base_aa_seq))

        # If we require more than 1 change
        if self.num_changes > 1:

            # Loop through how many changes we require
            for change_no in range(2, self.num_changes+1):

                # Loop through 
                for idx, seq in enumerate(self.generated_aa_seq):

                    # Only generate new sequences using existing generated sequences as a base, if it has undergone 1 change less than required
                    if sum(self.generated_aa_num_changes[idx]) == (change_no - 1):

                        # Invert generated_aa_num_changes list
                        change_mask = [not elem for elem in self.generated_aa_num_changes[idx]]

                        # Compute a new change vector by ANDing the change_mask and the change_aa_vector (to remove the amino acid(s) that have already changed)
                        new_aa_change_vec = [a and b for a, b in zip(change_mask, self.change_aa_vector[:])]
                        new_nt_change_vec, _ = ogu.sync_aa_change_to_nt_change(new_aa_change_vec, self.change_nt_vector[:])

                        # Generate new sequences
                        new_aa_seq, new_aa_num_changes = ogu.generate_aa_sequences(self.base_aa_seq, self.base_nt_seq, new_aa_change_vec, new_nt_change_vec, self.generated_aa_num_changes[idx])

                        # Get unique values by using the set() operator
                        unique_idx         = [new_aa_seq.index(x) for x in set(new_aa_seq)]
                        unique_seq         = [new_aa_seq[idx] for idx in unique_idx]
                        unique_num_changes = [new_aa_num_changes[idx] for idx in unique_idx]

                        # Append unique values to internal properties
                        self.generated_aa_seq         = self.generated_aa_seq + unique_seq
                        self.generated_aa_num_changes = self.generated_aa_num_changes + unique_num_changes

        return
    
    def generate_nt_sequences(self):

        # Generate the nucleotide sequences
        self.generated_nt_seq, self.generated_nt_seq_change_attempts = ogu.generate_nt_sequences(
            self.generated_aa_seq, self.generated_aa_num_changes, self.base_nt_seq, 
            self.change_nt_vector, self.fullyfree_vector, self.restriction_sites)

        return

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
        self._base_aa_seq = ogu.nt2aa( self._base_nt_seq )

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
        self._change_aa_vector, self._fullyfree_vector = ogu.sync_nt_change_to_aa_change(value, self._change_aa_vector)
    
    @property
    def change_aa_vector(self):
        return self._change_aa_vector
    
    @change_aa_vector.setter
    def change_aa_vector(self, value):
        # Assign value
        self._change_aa_vector = value
        # Sync values to the nucleotide vector
        self._change_nt_vector, self._fullyfree_vector = ogu.sync_aa_change_to_nt_change(value, self._change_nt_vector)
    
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
        