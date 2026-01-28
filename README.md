# SLOGGER

**Saturation Libraries, Golden-Gate Excision Ready**

SLOGGER is a deterministic, Dockerized pipeline for generating site-saturation mutagenesis libraries from protein or coding sequences and producing Golden Gate-ready DNA fragments. It is designed for reproducibility, auditability, and tight integration with codon optimization and cloning workflows.

SLOGGER generates exhaustive single-site amino-acid substitution libraries (including wild type), validates post-ligation ORFs, and outputs fragment FASTAs suitable for direct ordering or downstream assembly.

---

## Key Features

- Full site-saturation mutagenesis
  - Every amino acid position
  - All 20 amino acids (WT + 19 substitutions)
  - Deterministic ordering and outputs
- Codon optimization first
  - Uses `codonopt` as the baseline DNA generator
  - Same codon tables, constraints, and optimization modes
  - Codon selection at mutation sites respects forbidden motifs
- Golden Gate-ready outputs
  - User-defined fragment cut points
  - User-supplied adapters
  - Fragment FASTA is the primary output
- Strict post-ligation validation
  - Translation correctness
  - Length multiple of 3
  - No internal stop codons
  - GC bounds and homopolymer limits
  - Forbidden motif detection with context
- Deterministic and auditable
  - Fixed seeds
  - Stable variant ordering
  - Complete TSV report of all variants and failures
- Dockerized
  - No conda
  - Reproducible builds
  - Designed for HPC and workstation use

---

## What SLOGGER Produces

For each input gene, SLOGGER outputs:

- `baseline_optimized.fasta`  
  Codon-optimized baseline CDS (one per gene)
- `variants.fasta`  
  Full-length variant CDS sequences (adapter-free)
- `variants_protein.fasta`  
  Protein sequences corresponding to each variant
- `fragments_with_gg.fasta` (primary output)  
  Golden Gate-ready fragments with adapters
- `reassembled_orfs.fasta`  
  Adapter-free post-ligation ORFs used for validation
- `report.tsv`  
  Variant-level audit table (pass/fail, GC, motifs, etc.)
- `junction_hits.tsv`  
  Detailed forbidden-motif hit locations with sequence context

---

## Installation

### Clone with submodules

SLOGGER depends on `codonopt`, which is included as a git submodule.

```bash
git clone --recurse-submodules https://github.com/YOUR-USERNAME/slogger.git
cd slogger
```

If you already cloned without submodules:

```bash
git submodule update --init --recursive
```

---

## Docker (Recommended)

### Build

```bash
docker build -t slogger:0.1 .
```

### Run

```bash
docker run --rm \
  -v "$PWD:/work" \
  -w /work \
  slogger:0.1 \
  run config.yaml --outdir out
```

All outputs will be written to `out/`.

---

## Configuration

SLOGGER is driven by a single YAML config file.

### Minimal Example

```yaml
pipeline_version: "0.1"

inputs:
  fasta: inputs.fasta
  input_type: protein

codonopt:
  sequence: inputs.fasta
  codon_table: codon_table.xlsx
  codon_table_sheet: Ecoli
  out_subdir: codonopt_out
  gc_min: 0.35
  gc_max: 0.65
  max_homopolymer: 5
  seed: 1337

splits:
  global_cut_points: [300, 600]
  enforce_frame: true

golden_gate:
  fragments:
    - left_adapter: GGTCTC
      right_adapter: GAGACC
    - left_adapter: GGTCTC
      right_adapter: GAGACC
    - left_adapter: GGTCTC
      right_adapter: GAGACC

final_checks:
  require_no_internal_stops: true
  require_length_multiple_of_3: true

outputs:
  fragments_fasta: fragments_with_gg.fasta
  reassembled_fasta: reassembled_orfs.fasta
  report_tsv: report.tsv
```

---

## Mutation Model

- Single-site saturation only
- One position mutated at a time
- WT amino acid is included explicitly
- Codons at mutated sites are selected deterministically:
  1. Highest-frequency codon for the target amino acid
  2. If that introduces a forbidden motif, the next-best synonymous codon is tried
- No random sampling is performed

---

## Fragmentation and Adapters

- Fragment cut points are 0-based nucleotide indices
- Frame enforcement is optional
- Adapter sequences are trusted verbatim
- Validation is performed on the adapter-free reassembled ORF

---

## Determinism and Reproducibility

Given the same:
- input FASTA
- YAML config
- codon table
- codonopt commit

SLOGGER will always produce identical outputs.

---

## Development (Optional)

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest
```

---

## When to Use SLOGGER

SLOGGER is ideal when you want:
- Complete saturation mutagenesis libraries
- Codon-optimized, synthesis-ready designs
- Golden Gate cloning outputs
- Full traceability and validation
- Deterministic, repeatable builds
