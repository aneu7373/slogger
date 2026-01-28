# SLOGGER
**S**aturation **L**ibraries, **G**olden-Gate **E**xcision **R**eady

SLOGGER is a deterministic site-saturation mutagenesis design pipeline that:
- runs **codonopt** to generate a baseline codon-optimized CDS per gene
- generates a complete saturation library (WT + 19 substitutions) at every amino-acid position
- **always** splits variants into fragments using FROGGER-compatible cut points
- **always** adds Golden Gate adapters (output-only)
- validates constructs using FROGGER-compatible post-ligation ORF checks

## Repo setup (with codonopt as submodule)

From the repo root:

```bash
git init
git add .
git commit -m "Initial SLOGGER scaffold"

# Add codonopt as a submodule (same pattern as FROGGER)
git submodule add <YOUR_CODONOPT_GIT_URL> third_party/codonopt
git submodule update --init --recursive

git commit -m "Add codonopt submodule"
```

> Replace `<YOUR_CODONOPT_GIT_URL>` with the URL for your codonopt repository.

## Local install (dev)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -U pip

# install codonopt editable (from submodule)
pip install -e third_party/codonopt

# install slogger editable
pip install -e ".[dev]"

pytest -q
```

## Running SLOGGER

SLOGGER takes a FROGGER-compatible YAML config (see `FROGGER.schema.md` for the contract). At minimum you need:
- `inputs.fasta`
- `codonopt` block (sequence or batch_table, codon_table, etc.)
- `splits` cut points
- `golden_gate.fragments` adapters

Example:
```bash
slogger run config.yaml --outdir out/
```

Outputs (typical):
- `out/baseline_optimized.fasta`
- `out/variants.fasta`
- `out/variants_protein.fasta`
- `out/fragments_with_gg.fasta`
- `out/reassembled_orfs.fasta`
- `out/report.tsv`
- `out/junction_hits.tsv`

## Docker build & run

### Build (important: include submodules)
If you cloned the repo, ensure submodules are present:

```bash
git submodule update --init --recursive
```

Then build:

```bash
docker build -t slogger:0.1 .
```

### Run
Mount your config + inputs into the container:

```bash
docker run --rm \
  -v "$PWD:/work" \
  -w /work \
  slogger:0.1 \
  run config.yaml --outdir out
```

## Ensuring codonopt integration matches FROGGER

SLOGGER runs codonopt via:

```text
python -m codonopt ...
```

The YAML keys used to build the command come from `codonopt` in your config. If your codonopt CLI flags differ from what SLOGGER assumes, update `src/slogger/pipeline.py::run_codonopt()` to match the exact invocation used by FROGGER.
