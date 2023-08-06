# import statements
import python_codon_tables as pct
import oligo_generator_utility as ogu

class oligo_generator:

    def __init__(self, base_nt_seq='', num_changes=1):
        
        # Assign properties
        self.base_nt_seq = base_nt_seq
        self.num_changes = num_changes



    def generate_aa_sequences(self):

        # create the array of new amino acid sequences
        self.generated_aa_seq = ogu.generate_aa_sequences(self.base_aa_seq, self.change_aa_vector)

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
        self._change_nt_vector = [0] * len(self.base_nt_seq)
        self._change_aa_vector = [0] * len(self.base_aa_seq)
        self._fullyfree_vector = [0] * len(self.base_aa_seq)

    @property
    def base_aa_seq(self):
        return self._base_aa_seq
    
    @property
    def change_nt_vector(self):
        return self._change_nt_vector
    
    @property
    def change_aa_vector(self):
        return self._change_aa_vector
    
    @change_aa_vector.setter
    def change_aa_vector(self, value):
        # Assign value
        self._change_aa_vector = value
    
    @property
    def fullyfree_vector(self):
        return self._fullyfree_vector
        