from __future__ import annotations
import os
import glob
import subprocess
from typing import Any, Dict, List, Tuple

import yaml

from .io_fasta import read_fasta, write_fasta
from .codon_table_xlsx import load_codon_table_xlsx
from .mutate import generate_saturation_variants
from .fragments.frogger_backend import FroggerBackendV0_1


def load_config(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    if not isinstance(cfg, dict):
        raise ValueError("Config must be a YAML mapping.")
    return cfg


def _merge_cut_points(cfg: Dict[str, Any], gene_id: str) -> Tuple[List[int], bool]:
    splits = cfg.get("splits", {}) or {}
    global_cps = list(splits.get("global_cut_points", []) or [])
    per_gene = (splits.get("per_gene_cut_points", {}) or {})
    enforce = bool(splits.get("enforce_frame", True))
    if gene_id in per_gene and per_gene[gene_id] is not None:
        return list(per_gene[gene_id]), enforce
    return global_cps, enforce


def _collect_forbidden_motifs(cfg: Dict[str, Any]) -> List[str]:
    codonopt_cfg = cfg.get("codonopt", {}) or {}
    avoid = codonopt_cfg.get("avoid_motifs", []) or []
    final_forbidden = ((cfg.get("final_forbidden", {}) or {}).get("motifs", []) or [])
    out: List[str] = []
    out += [str(m).upper() for m in avoid if str(m).strip()]
    out += [str(m).upper() for m in final_forbidden if str(m).strip()]
    return out


def _find_codonopt_output_fasta(out_dir: str) -> str:
    # codonopt CLI default --out is "results" and historically many pipelines write optimized_sequences.fasta
    candidates = [
        os.path.join(out_dir, "optimized_sequences.fasta"),
        os.path.join(out_dir, "optimized_sequences.fa"),
        os.path.join(out_dir, "optimized.fasta"),
        os.path.join(out_dir, "results", "optimized_sequences.fasta"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c

    fastas = sorted(glob.glob(os.path.join(out_dir, "**", "*.fasta"), recursive=True))
    fastas += sorted(glob.glob(os.path.join(out_dir, "**", "*.fa"), recursive=True))
    fastas = [f for f in fastas if os.path.isfile(f)]

    uniq: List[str] = []
    for f in fastas:
        if f not in uniq:
            uniq.append(f)

    if len(uniq) == 1:
        return uniq[0]
    raise FileNotFoundError(
        f"Could not uniquely identify codonopt output FASTA in {out_dir}. Found: {uniq}"
    )


def _inject_codonopt_pythonpath(env: Dict[str, str], workdir: str) -> Dict[str, str]:
    """
    Ensure codonopt is importable even if it is NOT pip-installable.
    Expected submodule location: <repo_root>/third_party/codonopt
    """
    candidates = [
        os.path.join(workdir, "third_party", "codonopt"),
        os.path.join(os.path.dirname(workdir), "third_party", "codonopt"),
    ]
    codonopt_root = None
    for c in candidates:
        if os.path.isdir(os.path.join(c, "codonopt")):
            codonopt_root = c
            break

    if codonopt_root:
        existing = env.get("PYTHONPATH", "")
        env["PYTHONPATH"] = f"{codonopt_root}{os.pathsep}{existing}" if existing else codonopt_root

    return env


def run_codonopt(cfg: Dict[str, Any], workdir: str) -> str:
    """
    Run codonopt and return path to optimized CDS FASTA.

    This matches codonopt's CLI flags from codonopt/cli.py:
      --sequence / --batch-table
      --codon-table
      --codon-table-sheet
      --out
      --avoid-codons (pipe-separated string)
      --avoid-motifs (pipe-separated string)
      --gc-min/--gc-max
      --max-homopolymer
      --min-codon-fraction
      --seed
      --n
      --optimization-mode
      --max-tries-per-replicate
      --kleinbub-search-limit
      --backtrack-window
      --cryptkeeper-enable (+ related cryptkeeper flags)
    """
    codonopt_cfg = cfg.get("codonopt", {}) or {}

    # Where codonopt should write its outputs
    out_subdir = codonopt_cfg.get("out_subdir", "codonopt_out")
    out_dir = os.path.join(workdir, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    # Inputs: choose ONE
    seq = codonopt_cfg.get("sequence", None)
    batch = codonopt_cfg.get("batch_table", None)
    if bool(seq) == bool(batch):
        raise ValueError("codonopt: provide exactly one of 'sequence' or 'batch_table'.")

    codon_table = codonopt_cfg.get("codon_table", None)
    sheet = codonopt_cfg.get("codon_table_sheet", None)
    if not codon_table:
        raise ValueError("codonopt.codon_table is required.")

    # Codonopt entrypoint: module execution (requires PYTHONPATH to include third_party/codonopt)
    # If your codonopt repo uses a different executable entrypoint, change ONLY this.
    cmd: List[str] = ["python", "-m", "codonopt.main"]

    if seq:
        cmd += ["--sequence", str(seq)]
    if batch:
        cmd += ["--batch-table", str(batch)]

    cmd += ["--codon-table", str(codon_table)]
    if sheet:
        cmd += ["--codon-table-sheet", str(sheet)]

    # IMPORTANT: codonopt uses --out for output directory
    cmd += ["--out", out_dir]

    # Logging knobs
    if codonopt_cfg.get("log_file"):
        cmd += ["--log-file", str(codonopt_cfg["log_file"])]
    if bool(codonopt_cfg.get("verbose", False)):
        cmd += ["--verbose"]

    # Avoid lists become pipe-separated strings (per codonopt CLI)
    avoid_codons = codonopt_cfg.get("avoid_codons", []) or []
    avoid_motifs = codonopt_cfg.get("avoid_motifs", []) or []
    if avoid_codons:
        cmd += ["--avoid-codons", "|".join(map(str, avoid_codons))]
    if avoid_motifs:
        cmd += ["--avoid-motifs", "|".join(map(str, avoid_motifs))]

    # Basic constraints / knobs
    for k, flag in [
        ("gc_min", "--gc-min"),
        ("gc_max", "--gc-max"),
        ("max_homopolymer", "--max-homopolymer"),
        ("min_codon_fraction", "--min-codon-fraction"),
        ("seed", "--seed"),
        ("n", "--n"),
    ]:
        v = codonopt_cfg.get(k, None)
        if v is not None:
            cmd += [flag, str(v)]

    # Optimization behavior
    opt = codonopt_cfg.get("optimization", {}) or {}
    if opt.get("mode"):
        cmd += ["--optimization-mode", str(opt["mode"])]
    if opt.get("max_tries_per_replicate") is not None:
        cmd += ["--max-tries-per-replicate", str(opt["max_tries_per_replicate"])]
    if opt.get("kleinbub_search_limit") is not None:
        cmd += ["--kleinbub-search-limit", str(opt["kleinbub_search_limit"])]
    if opt.get("backtrack_window") is not None:
        cmd += ["--backtrack-window", str(opt["backtrack_window"])]

    # CryptKeeper flags (match codonopt CLI)
    ck = codonopt_cfg.get("cryptkeeper", {}) or {}
    if bool(ck.get("enable", False)):
        cmd += ["--cryptkeeper-enable"]
        if ck.get("rbs_score_cutoff") is not None:
            cmd += ["--cryptkeeper-rbs-score-cutoff", str(ck["rbs_score_cutoff"])]
        if ck.get("threads") is not None:
            cmd += ["--cryptkeeper-threads", str(ck["threads"])]
        if ck.get("fail_score") is not None:
            cmd += ["--cryptkeeper-fail-score", str(ck["fail_score"])]
        if ck.get("ignore_first_nt") is not None:
            cmd += ["--cryptkeeper-ignore-first-nt", str(ck["ignore_first_nt"])]
        if ck.get("ignore_last_nt") is not None:
            cmd += ["--cryptkeeper-ignore-last-nt", str(ck["ignore_last_nt"])]
        if ck.get("pool_factor") is not None:
            cmd += ["--cryptkeeper-pool-factor", str(ck["pool_factor"])]
        if ck.get("max_pool") is not None:
            cmd += ["--cryptkeeper-max-pool", str(ck["max_pool"])]

    env = dict(os.environ)
    env = _inject_codonopt_pythonpath(env, workdir)

    proc = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "codonopt failed.\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}\n"
            f"PYTHONPATH:\n{env.get('PYTHONPATH','')}\n"
        )

    return _find_codonopt_output_fasta(out_dir)


def run_pipeline(config_path: str, outdir: str) -> None:
    cfg = load_config(config_path)
    os.makedirs(outdir, exist_ok=True)

    inputs = cfg.get("inputs", {}) or {}
    fasta_in = inputs.get("fasta")
    if not fasta_in:
        raise ValueError("inputs.fasta is required.")

    # Run directory is the config directory
    workdir = os.path.dirname(os.path.abspath(config_path)) or "."

    # 1) codonopt baseline
    optimized_fasta = run_codonopt(cfg, workdir=workdir)
    optimized = read_fasta(optimized_fasta)

    # 2) load codon table used by codonopt (xlsx)
    codonopt_cfg = cfg.get("codonopt", {}) or {}
    table = load_codon_table_xlsx(codonopt_cfg["codon_table"], codonopt_cfg.get("codon_table_sheet"))

    forbidden = _collect_forbidden_motifs(cfg)

    backend = FroggerBackendV0_1()

    # outputs
    baseline_fasta = os.path.join(outdir, "baseline_optimized.fasta")
    variants_fasta = os.path.join(outdir, "variants.fasta")
    variants_prot_fasta = os.path.join(outdir, "variants_protein.fasta")

    fragments_fasta = os.path.join(
        outdir, (cfg.get("outputs", {}) or {}).get("fragments_fasta", "fragments_with_gg.fasta")
    )
    reassembled_fasta = os.path.join(
        outdir, (cfg.get("outputs", {}) or {}).get("reassembled_fasta", "reassembled_orfs.fasta")
    )
    report_tsv = os.path.join(
        outdir, (cfg.get("outputs", {}) or {}).get("report_tsv", "report.tsv")
    )
    hits_tsv = os.path.join(
        outdir, (cfg.get("outputs", {}) or {}).get("junction_hits_tsv", "junction_hits.tsv")
    )

    # write baseline
    write_fasta(baseline_fasta, [(r.id, r.seq) for r in optimized])

    # golden gate adapter config (list of dicts)
    gg_cfg = ((cfg.get("golden_gate", {}) or {}).get("fragments", []) or [])
    final_checks = cfg.get("final_checks", {}) or {}
    final_forbidden = (cfg.get("final_forbidden", {}) or {}).get("motifs", []) or []

    report_rows: List[Dict[str, str]] = []
    hit_rows: List[Dict[str, str]] = []

    variant_fasta_entries = []
    variant_prot_entries = []
    fragment_entries = []
    reassembled_entries = []

    # For each gene record from codonopt output
    for rec in optimized:
        gene_id = rec.id
        baseline_cds = rec.seq

        cut_points, enforce_frame = _merge_cut_points(cfg, gene_id)

        for var in generate_saturation_variants(
            gene_id=gene_id,
            baseline_cds=baseline_cds,
            forbidden_motifs=forbidden,
            codon_table=table,
            include_wt=True,
        ):
            variant_id = f"{gene_id}|{var.mutation_label}"

            # Split + adapters (FROGGER-compatible)
            frags_core = backend.split_sequence_by_cut_points(
                var.cds, cut_points=cut_points, enforce_frame=enforce_frame
            )
            frags_clone = backend.add_adapters(frags_core, gg_cfg=gg_cfg)

            # Reassemble core ORF
            core_orf = backend.reassemble_core_orf(frags_core)

            # Validate post-ligation ORF (core-only)
            vres = backend.check_construct(
                core_orf_seq=core_orf,
                codonopt_cfg=codonopt_cfg,
                final_forbidden_motifs=final_forbidden,
                require_no_internal_stops=bool(final_checks.get("require_no_internal_stops", True)),
                require_len_multiple_of_3=bool(final_checks.get("require_length_multiple_of_3", True)),
            )

            # FASTA outputs
            variant_fasta_entries.append((variant_id, var.cds))
            variant_prot_entries.append((variant_id, var.protein))
            reassembled_entries.append((variant_id, core_orf))
            for f in frags_clone:
                header = f"{variant_id}|frag{f.frag_index:02d}"
                fragment_entries.append((header, f.cloning_seq))

            # Report
            report_rows.append(
                {
                    "gene_id": gene_id,
                    "variant_id": variant_id,
                    "position_aa": str(var.pos1),
                    "wt_aa": var.wt_aa,
                    "target_aa": var.target_aa,
                    "mutation": var.mutation_label,
                    "codon_used": var.codon_used,
                    "passes_all": str(vres.passes_all),
                    "fail_reasons": ";".join(vres.fail_reasons),
                    "gc_fraction": str(vres.gc_fraction),
                    "max_homopolymer": str(vres.max_homopolymer),
                    "aa_len": str(vres.aa_len),
                    "n_frags": str(len(frags_core)),
                }
            )

            for h in vres.hits:
                hit_rows.append(
                    {
                        "variant_id": variant_id,
                        "motif": h["motif"],
                        "start": str(h["start"]),
                        "end": str(h["end"]),
                        "context": h["context"],
                    }
                )

    write_fasta(variants_fasta, variant_fasta_entries)
    write_fasta(variants_prot_fasta, variant_prot_entries)
    write_fasta(fragments_fasta, fragment_entries)
    write_fasta(reassembled_fasta, reassembled_entries)

    # TSVs
    def write_tsv(path: str, rows: List[Dict[str, str]], header: List[str]) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(r.get(h, "") for h in header) + "\n")

    write_tsv(
        report_tsv,
        report_rows,
        header=[
            "gene_id",
            "variant_id",
            "position_aa",
            "wt_aa",
            "target_aa",
            "mutation",
            "codon_used",
            "passes_all",
            "fail_reasons",
            "gc_fraction",
            "max_homopolymer",
            "aa_len",
            "n_frags",
        ],
    )
    write_tsv(hits_tsv, hit_rows, header=["variant_id", "motif", "start", "end", "context"])
