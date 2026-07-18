from pathlib import Path
from collections.abc import Iterable
from itertools import combinations
from typing import Any

import oligo_generator.utilities.helper_functions as ogu


DNA_BASES = frozenset("ACGT")


class oligo_generator:
    """Generate amino-acid and nucleotide oligo libraries.

    Public position arguments are 1-based so they match biological convention
    and the web UI. Internally, positions are converted to Python indices.
    """

    def __init__(self, base_nt_seq: str = '', num_changes: int = 1) -> None:
        """Create a generator from a base nucleotide sequence.

        Args:
            base_nt_seq: DNA sequence using A/C/G/T. Empty is allowed for
                construction, but generation requires at least one codon.
            num_changes: Maximum number of amino-acid substitutions to include.

        Raises:
            TypeError: If argument types are invalid.
            ValueError: If `base_nt_seq` contains invalid bases or has a length
                that is not divisible by three.
        """

        self.base_nt_seq = base_nt_seq
        self.num_changes = num_changes

        self.generated_aa_seq: list[str] = []
        self.generated_nt_seq: list[str] = []
        self.generated_aa_changes: list[list[bool]] = []
        self.generated_nt_seq_change_attempts: list[list[int]] = []
        self.restriction_sites: list[str] = []
        self.restricted_aa_sequences: list[str] = []

    def generate_aa_sequences(self, s_io: Any = None) -> None:
        """Generate amino-acid variants up to `num_changes` substitutions.

        Args:
            s_io: Optional SocketIO-like object with an `emit` method for
                progress updates.

        Raises:
            ValueError: If the base sequence is empty or `num_changes` exceeds
                the amino-acid sequence length.
            TypeError: If helper output types are not as expected.
        """

        self._validate_generation_ready()

        mutable_positions = [
            idx for idx, can_change in enumerate(self.change_aa_vector)
            if can_change
        ]

        for change in range(1, self.num_changes + 1):
            print('Generating sequences for ', change, ' number of changes.')

            change_combos = []
            for changed_positions in combinations(mutable_positions, change):
                combo = [False] * len(self.base_aa_seq)
                for idx in changed_positions:
                    combo[idx] = True
                change_combos.append(combo)

            generated_aa_seq, generated_aa_changes = ogu.generate_aa_sequences(
                self.base_aa_seq,
                change_combos,
                self.restricted_aa_sequences,
            )
            self._validate_aa_generation_output(generated_aa_seq,
                                                generated_aa_changes)

            self.generated_aa_seq += generated_aa_seq
            self.generated_aa_changes += generated_aa_changes

            if s_io is not None:
                progress = (change / self.num_changes) / 2 * 100
                s_io.emit('update_progress', {
                    'progress': progress,
                    'current_state': change,
                    'total': self.num_changes,
                })

    def generate_nt_sequences(self, s_io: Any = None) -> None:
        """Generate nucleotide sequences for generated amino-acid variants.

        Args:
            s_io: Optional SocketIO-like object with an `emit` method for
                progress updates.

        Raises:
            ValueError: If amino-acid variants have not been generated first.
            TypeError: If helper output types are not as expected.
        """

        if not self.generated_aa_seq or not self.generated_aa_changes:
            raise ValueError('Call generate_aa_sequences() before generating '
                             'nucleotide sequences.')

        generated_nt, change_attempts, kept_aa_indices = ogu.generate_nt_sequences(
            self.generated_aa_seq,
            self.generated_aa_changes,
            self.base_nt_seq,
            self.change_nt_vector,
            self.fullyfree_vector,
            self.restriction_sites,
            s_io=s_io,
        )
        self._validate_nt_generation_output(generated_nt, change_attempts,
                                            kept_aa_indices)

        self.generated_nt_seq = generated_nt
        self.generated_nt_seq_change_attempts = change_attempts
        self.generated_aa_seq = [self.generated_aa_seq[idx]
                                 for idx in kept_aa_indices]
        self.generated_aa_changes = [self.generated_aa_changes[idx]
                                     for idx in kept_aa_indices]

    def set_fixed_aa(self, value: str) -> None:
        """Fix every occurrence of an amino-acid code in the base sequence.

        Args:
            value: Single-letter amino-acid code to keep unchanged.
        """

        value = self._validate_amino_acid_code(value)
        for idx, aa in enumerate(self.base_aa_seq):
            if aa == value:
                self._set_aa_index(idx, False)

    def set_aa_pos(self, position: int, value: bool) -> None:
        """Set whether a 1-based amino-acid position may change.

        Args:
            position: 1-based amino-acid position.
            value: True if the amino acid can change, False to fix it.
        """

        value = self._validate_bool(value, 'value')
        idx = self._position_to_index(position, len(self.base_aa_seq),
                                      'Amino acid')
        self._set_aa_index(idx, value)

    def set_nt_pos(self, position: int, value: bool) -> None:
        """Set whether a 1-based nucleotide position may change.

        Args:
            position: 1-based nucleotide position.
            value: True if the nucleotide can change, False to fix it.
        """

        value = self._validate_bool(value, 'value')
        idx = self._position_to_index(position, len(self.base_nt_seq),
                                      'Nucleotide')
        self._set_nt_index(idx, value)

    def export_aa_sequences(self, fn: str | Path) -> None:
        """Write generated amino-acid sequences to a newline-delimited file."""

        path = self._validate_export_path(fn)
        if self.generated_aa_seq:
            path.write_text('\n'.join(self.generated_aa_seq))

    def export_nt_sequences(self, fn: str | Path) -> None:
        """Write generated nucleotide sequences to a newline-delimited file."""

        path = self._validate_export_path(fn)
        if self.generated_nt_seq:
            path.write_text('\n'.join(self.generated_nt_seq))

    @property
    def base_nt_seq(self) -> str:
        """Base DNA sequence used to derive amino-acid variants."""

        return self._base_nt_seq

    @base_nt_seq.setter
    def base_nt_seq(self, value: str) -> None:
        value = self._validate_dna_sequence(value, 'base_nt_seq')
        self._base_nt_seq = value
        self._base_aa_seq = ogu.nt2aa(self._base_nt_seq)
        self._change_nt_vector = [True] * len(self.base_nt_seq)
        self._change_aa_vector = [True] * len(self.base_aa_seq)
        self._fullyfree_vector = [True] * len(self.base_aa_seq)

    @property
    def base_aa_seq(self) -> str:
        """Amino-acid translation of `base_nt_seq`."""

        return self._base_aa_seq

    @property
    def num_changes(self) -> int:
        """Maximum number of amino-acid substitutions to generate."""

        return self._num_changes

    @num_changes.setter
    def num_changes(self, value: int) -> None:
        if isinstance(value, bool) or not isinstance(value, int):
            raise TypeError('num_changes must be an integer.')
        if value < 1:
            raise ValueError('num_changes must be at least 1.')
        self._num_changes = value

    @property
    def change_nt_vector(self) -> list[bool]:
        """Per-nucleotide mutability vector."""

        return self._change_nt_vector

    @change_nt_vector.setter
    def change_nt_vector(self, value: list[bool]) -> None:
        value = self._validate_bool_vector(value, len(self.base_nt_seq),
                                           'change_nt_vector')
        self._change_nt_vector = value
        self._change_aa_vector, self._fullyfree_vector = (
            ogu.sync_nt_change_to_aa_change(value, self._change_aa_vector)
        )
        self._validate_internal_vectors()

    @property
    def change_aa_vector(self) -> list[bool]:
        """Per-amino-acid mutability vector."""

        return self._change_aa_vector

    @change_aa_vector.setter
    def change_aa_vector(self, value: list[bool]) -> None:
        value = self._validate_bool_vector(value, len(self.base_aa_seq),
                                           'change_aa_vector')
        self._change_aa_vector = value
        self._change_nt_vector, self._fullyfree_vector = (
            ogu.sync_aa_change_to_nt_change(value, self._change_nt_vector)
        )
        self._validate_internal_vectors()

    @property
    def fullyfree_vector(self) -> list[bool]:
        """Per-amino-acid flag indicating all three codon bases are mutable."""

        return self._fullyfree_vector

    @fullyfree_vector.setter
    def fullyfree_vector(self, value: list[bool]) -> None:
        self._fullyfree_vector = self._validate_bool_vector(
            value, len(self.base_aa_seq), 'fullyfree_vector'
        )

    @property
    def restriction_sites(self) -> list[str]:
        """Nucleotide motifs that generated DNA sequences must exclude."""

        return self._restriction_sites

    @restriction_sites.setter
    def restriction_sites(self, value: Iterable[str]) -> None:
        self._restriction_sites = [
            self._validate_dna_motif(site, 'restriction site')
            for site in self._validate_sequence_collection(value,
                                                           'restriction_sites')
            if site
        ]

    @property
    def restricted_aa_sequences(self) -> list[str]:
        """Amino-acid motifs that generated protein sequences must exclude."""

        return self._restricted_aa_sequences

    @restricted_aa_sequences.setter
    def restricted_aa_sequences(self, value: Iterable[str]) -> None:
        sequences = self._validate_sequence_collection(
            value, 'restricted_aa_sequences'
        )
        self._restricted_aa_sequences = [seq.upper() for seq in sequences if seq]

    def _position_to_index(self, position: int, sequence_length: int,
                           sequence_name: str) -> int:
        if isinstance(position, bool) or not isinstance(position, int):
            raise TypeError(f'{sequence_name} position must be an integer.')

        idx = position - 1
        if idx < 0 or idx >= sequence_length:
            raise IndexError(
                f'{sequence_name} position {position} is out of range. '
                f'Use a 1-based position from 1 to {sequence_length}.'
            )

        return idx

    def _set_aa_index(self, idx: int, value: bool) -> None:
        aa_change_vec = self.change_aa_vector[:]
        aa_change_vec[idx] = value
        self.change_aa_vector = aa_change_vec

    def _set_nt_index(self, idx: int, value: bool) -> None:
        nt_change_vec = self.change_nt_vector[:]
        nt_change_vec[idx] = value
        self.change_nt_vector = nt_change_vec

    def _validate_generation_ready(self) -> None:
        if not self.base_aa_seq:
            raise ValueError('base_nt_seq must contain at least one codon.')
        if self.num_changes > len(self.base_aa_seq):
            raise ValueError('num_changes cannot exceed the amino-acid '
                             'sequence length.')

    def _validate_internal_vectors(self) -> None:
        self._validate_bool_vector(self._change_nt_vector,
                                   len(self.base_nt_seq),
                                   'change_nt_vector')
        self._validate_bool_vector(self._change_aa_vector,
                                   len(self.base_aa_seq),
                                   'change_aa_vector')
        self._validate_bool_vector(self._fullyfree_vector,
                                   len(self.base_aa_seq),
                                   'fullyfree_vector')

    @staticmethod
    def _validate_bool(value: bool, name: str) -> bool:
        if not isinstance(value, bool):
            raise TypeError(f'{name} must be a boolean.')
        return value

    @staticmethod
    def _validate_bool_vector(value: list[bool], expected_length: int,
                              name: str) -> list[bool]:
        if not isinstance(value, list):
            raise TypeError(f'{name} must be a list of booleans.')
        if len(value) != expected_length:
            raise ValueError(f'{name} must contain {expected_length} values.')
        if any(not isinstance(item, bool) for item in value):
            raise TypeError(f'{name} must contain only booleans.')
        return value[:]

    @staticmethod
    def _validate_dna_sequence(value: str, name: str) -> str:
        if not isinstance(value, str):
            raise TypeError(f'{name} must be a string.')

        value = value.strip().upper()
        if len(value) % 3 != 0:
            raise ValueError(f'{name} length must be divisible by three.')
        invalid_bases = sorted(set(value) - DNA_BASES)
        if invalid_bases:
            raise ValueError(f'{name} contains invalid DNA bases: '
                             f'{", ".join(invalid_bases)}.')
        return value

    @staticmethod
    def _validate_dna_motif(value: str, name: str) -> str:
        if not isinstance(value, str):
            raise TypeError(f'{name} must be a string.')

        value = value.strip().upper()
        invalid_bases = sorted(set(value) - DNA_BASES)
        if invalid_bases:
            raise ValueError(f'{name} contains invalid DNA bases: '
                             f'{", ".join(invalid_bases)}.')
        return value

    @staticmethod
    def _validate_amino_acid_code(value: str) -> str:
        if not isinstance(value, str):
            raise TypeError('Amino-acid code must be a string.')

        value = value.strip().upper()
        if len(value) != 1 or not value.isalpha():
            raise ValueError('Amino-acid code must be a single letter.')
        return value

    @staticmethod
    def _validate_sequence_collection(value: Iterable[str], name: str) -> list[str]:
        if isinstance(value, str) or not isinstance(value, Iterable):
            raise TypeError(f'{name} must be an iterable of strings.')

        sequences = []
        for item in value:
            if not isinstance(item, str):
                raise TypeError(f'{name} must contain only strings.')
            sequences.append(item.strip().upper())
        return sequences

    @staticmethod
    def _validate_export_path(fn: str | Path) -> Path:
        if not isinstance(fn, (str, Path)):
            raise TypeError('Export path must be a string or pathlib.Path.')
        return Path(fn)

    def _validate_aa_generation_output(self, sequences: list[str],
                                       changes: list[list[bool]]) -> None:
        if not isinstance(sequences, list) or not all(
            isinstance(seq, str) and len(seq) == len(self.base_aa_seq)
            for seq in sequences
        ):
            raise TypeError('Generated amino-acid sequences must be strings '
                            'matching the base amino-acid length.')

        if not isinstance(changes, list) or len(changes) != len(sequences):
            raise TypeError('Generated amino-acid change masks must align with '
                            'generated amino-acid sequences.')

        for change_mask in changes:
            self._validate_bool_vector(change_mask, len(self.base_aa_seq),
                                       'generated_aa_changes')

    def _validate_nt_generation_output(self, sequences: list[str],
                                       attempts: list[list[int]],
                                       kept_indices: list[int]) -> None:
        if not isinstance(sequences, list) or not all(
            isinstance(seq, str) and len(seq) == len(self.base_nt_seq)
            for seq in sequences
        ):
            raise TypeError('Generated nucleotide sequences must be strings '
                            'matching the base nucleotide length.')

        if not isinstance(attempts, list) or len(attempts) != len(sequences):
            raise TypeError('Codon-attempt metadata must align with generated '
                            'nucleotide sequences.')

        for attempt_vector in attempts:
            if not isinstance(attempt_vector, list) or not all(
                not isinstance(item, bool) and isinstance(item, int) and item >= 0 for item in attempt_vector
            ):
                raise TypeError('Codon-attempt metadata must contain lists of '
                                'non-negative integers.')

        if not isinstance(kept_indices, list) or len(kept_indices) != len(sequences):
            raise TypeError('Kept amino-acid indices must align with generated '
                            'nucleotide sequences.')

        max_index = len(self.generated_aa_seq) - 1
        if any(
            isinstance(idx, bool) or not isinstance(idx, int) or
            idx < 0 or idx > max_index
            for idx in kept_indices
        ):
            raise ValueError('Kept amino-acid indices are out of range.')
