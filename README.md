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
- `pyproject.toml` and `uv.lock`: uv project metadata and pinned dependency set.

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
- `uv` for dependency management and command execution.

The project dependencies are declared in `pyproject.toml` and pinned in `uv.lock`. Core runtime packages include Flask, Flask-SocketIO, gevent, SymPy, and python-codon-tables. Development tools such as pytest and flake8 are in the `dev` dependency group.

## Setup

Install `uv` if it is not already available:

```bash
python -m pip install uv
```

Create or update the local virtual environment from `uv.lock`:

```bash
uv sync --all-groups
```

For runtime dependencies only, omit the development group:

```bash
uv sync --no-dev
```

## Run The Web App

Start the Flask-SocketIO app:

```bash
uv run python startup.py
```

For local development you can also run the app module directly:

```bash
uv run python oligo_generator/webapp/oligo_app.py
```

The app listens on `0.0.0.0:8000` when run from `oligo_app.py`.

## Use The Python API

```python
import oligo_generator.models.generator as og

base_seq = "AGAAGCTGCATT"
o = og.oligo_generator(base_nt_seq=base_seq, num_changes=2)

o.set_aa_pos(3, False)      # Fix amino-acid position 3, using 1-based indexing.
o.set_nt_pos(5, False)      # Fix nucleotide position 5, using 1-based indexing.
o.set_fixed_aa("C")         # Keep every C in the base amino-acid sequence fixed.
o.restriction_sites = ["GGATCC", "GGTCTC"]
o.restricted_aa_sequences = ["RS"]

o.generate_aa_sequences()
o.generate_nt_sequences()

print(o.base_aa_seq)
print(o.generated_aa_seq)
print(o.generated_nt_seq)
```

## Position Indexing

All public position inputs are 1-based. This applies to both the web form and the Python API. For example, `set_aa_pos(1, False)` fixes the first amino acid and `set_nt_pos(1, False)` fixes the first nucleotide. Internally the implementation converts these values to Python list indices.

## Tests And Checks

Run the test suite with:

```bash
uv run pytest -q
```

Run the same lint checks used by CI:

```bash
uv run flake8 ./oligo_generator/ --count --select=E9,F63,F7,F82 --show-source --statistics
uv run flake8 ./oligo_generator/ --count --exit-zero --max-complexity=10 --max-line-length=127 --statistics
```

The suite currently has 16 tests. Coverage includes amino-acid generation counts, nucleotide generation, fixed amino-acid and nucleotide positions, fixed amino-acid codes, restricted amino-acid sequences, longer sequences with multiple fixed positions, higher change counts, and regressions for amino-acid/nucleotide alignment.

## Dependency Locking And Deployment

Use `uv lock` after changing dependencies in `pyproject.toml`:

```bash
uv lock
```

The Azure deployment workflow exports runtime dependencies to `requirements.txt` with:

```bash
uv export --no-dev --no-emit-project --no-hashes --format requirements.txt --output-file requirements.txt
```

`poetry.lock` is no longer used; `uv.lock` is the source of pinned dependency versions.

## Notes And Limitations

- Input nucleotide sequences should be uppercase DNA sequences whose length is divisible by three.
- The generator currently uses the `h_sapiens_9606` codon table.
- Generated substitutions intentionally exclude `*`, `C`, and `M`.
- Public position inputs in both the Python API and web form are one-based.
- Very large libraries can grow quickly because generated sequence count scales combinatorially with sequence length and number of changes.
