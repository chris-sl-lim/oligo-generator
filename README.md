# Oligo Generator

A Flask application and Python library for generating oligo libraries from a base nucleotide sequence. The generator translates the input nucleotide sequence to amino acids, creates amino-acid variants with a requested number of changes, then selects compatible nucleotide codons while respecting fixed positions and restricted sequence constraints.

## What Is In This Repository

- `oligo_generator/models/generator.py`: main `oligo_generator` class used by scripts and the web app.
- `oligo_generator/utilities/helper_functions.py`: nucleotide/amino-acid translation, variant generation, codon filtering, and restriction-site handling.
- `oligo_generator/webapp/oligo_app.py`: Flask and Flask-SocketIO web application with `/`, `/checkInputs`, and `/generateSequences` routes.
- `oligo_generator/webapp/templates/index.html`: browser UI for entering a base sequence, fixed positions, fixed amino-acid codes, restricted amino-acid sequences, and restriction sites.
- `RestrictionSites.csv`: restriction enzyme site data used by the example workflow.
- `example_oligo.py`: command-line example of using the generator class directly.
- `startup.py`: app entry point for hosting environments.
- `tests/`: pytest coverage for amino-acid and nucleotide generation behavior.
- `pyproject.toml` and `poetry.lock`: Poetry project metadata and pinned dependency set.

## Features

- Translates nucleotide sequences to amino-acid sequences using the human codon table from `python-codon-tables`.
- Generates amino-acid variant libraries for one or more amino-acid changes.
- Excludes stop codons, cysteine, methionine, and the original amino acid from generated substitutions.
- Supports fixed amino-acid positions, fixed nucleotide positions, and fixed amino-acid codes.
- Supports restricted amino-acid subsequences and restricted nucleotide sequences such as restriction sites.
- Selects codons by codon usage order while preserving nucleotide-position constraints.
- Provides a Flask UI with progress updates via SocketIO.

## Requirements

- Python 3.10 or newer.
- Poetry for dependency management.

The project dependencies are declared in `pyproject.toml`. Core runtime packages include Flask, Flask-SocketIO, gevent, SymPy, and python-codon-tables.

## Setup

Install dependencies with Poetry:

```bash
poetry install --with dev
```

If you only need runtime dependencies, omit `--with dev`.

## Run The Web App

Start the Flask-SocketIO app:

```bash
poetry run python startup.py
```

For local development you can also run the app module directly:

```bash
poetry run python oligo_generator/webapp/oligo_app.py
```

The app listens on `0.0.0.0:8000` when run from `oligo_app.py`.

## Use The Python API

```python
import oligo_generator.models.generator as og

base_seq = "AGAAGCTGCATT"
o = og.oligo_generator(base_nt_seq=base_seq, num_changes=2)

o.set_aa_pos(2, False)      # Fix amino-acid position, using zero-based indexing.
o.set_nt_pos(4, False)      # Fix nucleotide position, using zero-based indexing.
o.set_fixed_aa("C")         # Keep every C in the base amino-acid sequence fixed.
o.restriction_sites = ["GGATCC", "GGTCTC"]
o.restricted_aa_sequences = ["RS"]

o.generate_aa_sequences()
o.generate_nt_sequences()

print(o.base_aa_seq)
print(o.generated_aa_seq)
print(o.generated_nt_seq)
```

The web UI accepts biological 1-based positions and converts them internally to the generator's zero-based indexing.

## Tests

Run the test suite with:

```bash
poetry run pytest -q
```

The tests cover amino-acid generation, nucleotide generation, fixed positions, fixed amino-acid codes, restricted amino-acid sequences, and regressions for amino-acid/nucleotide alignment.

## Notes And Limitations

- Input nucleotide sequences should be uppercase DNA sequences whose length is divisible by three.
- The generator currently uses the `h_sapiens_9606` codon table.
- Generated substitutions intentionally exclude `*`, `C`, and `M`.
- Direct Python API positions are zero-based. Web-form positions are one-based.
