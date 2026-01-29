from __future__ import annotations

import glob
import os
import subprocess
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import yaml

from .io_fasta import read_fasta, write_fasta
from .codon_table_xlsx import load_codon_table_xlsx
from .mutate import (
    generate_explicit_variants,
    generate_saturation_variants,
    translate_cds,
)
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

    Matches codonopt CLI from your codonopt/cli.py:
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

    out_subdir = codonopt_cfg.get("out_subdir", "codonopt_out")
    out_dir = os.path.join(workdir, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    seq = codonopt_cfg.get("sequence", None)
    batch = codonopt_cfg.get("batch_table", None)
    if bool(seq) == bool(batch):
        raise ValueError("codonopt: provide exactly one of 'sequence' or 'batch_table'.")

    codon_table = codonopt_cfg.get("codon_table", None)
    sheet = codonopt_cfg.get("codon_table_sheet", None)
    if not codon_table:
        raise ValueError("codonopt.codon_table is required.")

    cmd: List[str] = ["python", "-m", "codonopt.main"]

    if seq:
        cmd += ["--sequence", str(seq)]
    if batch:
        cmd += ["--batch-table", str(batch)]

    cmd += ["--codon-table", str(codon_table)]
    if sheet:
        cmd += ["--codon-table-sheet", str(sheet)]

    cmd += ["--out", out_dir]

    if codonopt_cfg.get("log_file"):
        cmd += ["--log-file", str(codonopt_cfg["log_file"])]
    if bool(codonopt_cfg.get("verbose", False)):
        cmd += ["--verbose"]

    avoid_codons = codonopt_cfg.get("avoid_codons", []) or []
    avoid_motifs = codonopt_cfg.get("avoid_motifs", []) or []
    if avoid_codons:
        cmd += ["--avoid-codons", "|".join(map(str, avoid_codons))]
    if avoid_motifs:
        cmd += ["--avoid-motifs", "|".join(map(str, avoid_motifs))]

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

    opt = codonopt_cfg.get("optimization", {}) or {}
    if opt.get("mode"):
        cmd += ["--optimization-mode", str(opt["mode"])]
    if opt.get("max_tries_per_replicate") is not None:
        cmd += ["--max-tries-per-replicate", str(opt["max_tries_per_replicate"])]
    if opt.get("kleinbub_search_limit") is not None:
        cmd += ["--kleinbub-search-limit", str(opt["kleinbub_search_limit"])]
    if opt.get("backtrack_window") is not None:
        cmd += ["--backtrack-window", str(opt["backtrack_window"])]

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


MutationSpec = Tuple[str, Union[List[int], List[str]]]
# ("saturation_positions", [1,2,3]) OR ("explicit", ["T24V","A100G"])


def _parse_mutation_spec(cfg: Dict[str, Any], baseline_protein: str) -> MutationSpec:
    """
    mutation:
      aa_range: [start, end]     # inclusive, 1-based
    OR
      aa_positions: [1, 5, 99]
    OR
      explicit: ["T24V", "A100G"]

    If mutation block missing -> default saturation across ALL positions.

    Returns:
      ("saturation_positions", positions_to_mutate) OR ("explicit", explicit_list)
    """
    aa_len = len(baseline_protein)
    mut = cfg.get("mutation", None)
    if mut is None:
        return ("saturation_positions", list(range(1, aa_len + 1)))

    if not isinstance(mut, dict):
        raise ValueError("mutation must be a mapping if provided.")

    aa_range = mut.get("aa_range", None)
    aa_positions = mut.get("aa_positions", None)
    explicit = mut.get("explicit", None)

    provided = [x is not None for x in (aa_range, aa_positions, explicit)]
    if sum(provided) > 1:
        raise ValueError("Provide only one of mutation.aa_range, mutation.aa_positions, or mutation.explicit.")

    if explicit is not None:
        if not isinstance(explicit, list) or not all(isinstance(x, str) for x in explicit):
            raise ValueError("mutation.explicit must be a list of strings like ['T24V', 'A100G'].")
        cleaned = [s.strip() for s in explicit if s.strip()]
        if not cleaned:
            raise ValueError("mutation.explicit was provided but empty.")
        return ("explicit", cleaned)

    if aa_range is not None:
        if (
            not isinstance(aa_range, list)
            or len(aa_range) != 2
            or not all(isinstance(x, int) for x in aa_range)
        ):
            raise ValueError("mutation.aa_range must be a list of two integers: [start, end].")
        start, end = aa_range
        if start < 1 or end < 1 or end < start:
            raise ValueError("mutation.aa_range must satisfy 1 <= start <= end.")
        if end > aa_len:
            raise ValueError(f"mutation.aa_range end={end} exceeds protein length {aa_len}.")
        return ("saturation_positions", list(range(start, end + 1)))

    if aa_positions is not None:
        if not isinstance(aa_positions, list) or not all(isinstance(x, int) for x in aa_positions):
            raise ValueError("mutation.aa_positions must be a list of 1-based integer positions.")
        if any(p < 1 or p > aa_len for p in aa_positions):
            bad = [p for p in aa_positions if p < 1 or p > aa_len]
            raise ValueError(f"mutation.aa_positions contains out-of-range positions for aa_len={aa_len}: {bad}")
        # de-dup preserve order
        seen = set()
        uniq = []
        for p in aa_positions:
            if p not in seen:
                uniq.append(p)
                seen.add(p)
        return ("saturation_positions", uniq)

    # mutation block exists but contains none of the recognized keys
    return ("saturation_positions", list(range(1, aa_len + 1)))


def run_pipeline(config_path: str, outdir: str) -> None:
    cfg = load_config(config_path)
    os.makedirs(outdir, exist_ok=True)

    workdir = os.path.dirname(os.path.abspath(config_path)) or "."

    # 1) codonopt baseline
    optimized_fasta = run_codonopt(cfg, workdir=workdir)
    optimized = read_fasta(optimized_fasta)

    # 2) codon table used for mutation-site codon selection
    codonopt_cfg = cfg.get("codonopt", {}) or {}
    table = load_codon_table_xlsx(codonopt_cfg["codon_table"], codonopt_cfg.get("codon_table_sheet"))

    forbidden = _collect_forbidden_motifs(cfg)
    backend = FroggerBackendV0_1()

    # Common outputs always written
    baseline_fasta = os.path.join(outdir, "baseline_optimized.fasta")
    variants_fasta = os.path.join(outdir, "variants.fasta")
    variants_prot_fasta = os.path.join(outdir, "variants_protein.fasta")

    # Mode selection
    library_mode = (cfg.get("library_mode", {}) or {}).get("type", "full_orf_variants")
    max_frags_cfg = int((cfg.get("library_mode", {}) or {}).get("max_fragments", 8))

    # Golden Gate + checks config
    gg_cfg = ((cfg.get("golden_gate", {}) or {}).get("fragments", []) or [])
    final_checks = cfg.get("final_checks", {}) or {}
    final_forbidden = (cfg.get("final_forbidden", {}) or {}).get("motifs", []) or []

    # Output names
    out_cfg = cfg.get("outputs", {}) or {}
    fragments_fasta = os.path.join(outdir, out_cfg.get("fragments_fasta", "fragments_with_gg.fasta"))
    reassembled_fasta = os.path.join(outdir, out_cfg.get("reassembled_fasta", "reassembled_orfs.fasta"))
    report_tsv = os.path.join(outdir, out_cfg.get("report_tsv", "report.tsv"))
    hits_tsv = os.path.join(outdir, out_cfg.get("junction_hits_tsv", "junction_hits.tsv"))

    # New outputs for fragment pooling mode
    wt_frags_fasta = os.path.join(outdir, out_cfg.get("wt_fragments_fasta", "wt_fragments_with_gg.fasta"))
    mut_frags_fasta = os.path.join(outdir, out_cfg.get("mut_fragments_fasta", "mut_fragments_with_gg.fasta"))
    pools_dir = os.path.join(outdir, out_cfg.get("pools_dir", "pools"))
    assembly_plan_tsv = os.path.join(outdir, out_cfg.get("assembly_plan_tsv", "assembly_plan.tsv"))

    # Write baseline
    write_fasta(baseline_fasta, [(r.id, r.seq) for r in optimized])

    # Collect TSV rows
    report_rows: List[Dict[str, str]] = []
    hit_rows: List[Dict[str, str]] = []
    assembly_rows: List[Dict[str, str]] = []

    # Collect FASTA entries
    variant_fasta_entries: List[Tuple[str, str]] = []
    variant_prot_entries: List[Tuple[str, str]] = []
    reassembled_entries: List[Tuple[str, str]] = []

    # Fragment pooling outputs
    wt_fragment_entries: List[Tuple[str, str]] = []
    mut_fragment_entries: List[Tuple[str, str]] = []
    pools: Dict[int, List[Tuple[str, str]]] = {}  # frag_index -> fasta entries

    # Classic fragments output (full_orf_variants mode)
    classic_fragment_entries: List[Tuple[str, str]] = []

    def write_tsv(path: str, rows: List[Dict[str, str]], header: List[str]) -> None:
        with open(path, "w", encoding="utf-8") as f:
            f.write("\t".join(header) + "\n")
            for r in rows:
                f.write("\t".join(r.get(h, "") for h in header) + "\n")

    # --- MAIN LOOP ---
    for rec in optimized:
        gene_id = rec.id
        baseline_cds = rec.seq
        baseline_prot = translate_cds(baseline_cds)

        mut_mode, mut_payload = _parse_mutation_spec(cfg, baseline_prot)

        cut_points, enforce_frame = _merge_cut_points(cfg, gene_id)

        # WT fragments (core + cloning) from baseline
        wt_frags_core = backend.split_sequence_by_cut_points(
            baseline_cds, cut_points=cut_points, enforce_frame=enforce_frame
        )
        if len(wt_frags_core) > max_frags_cfg:
            raise ValueError(
                f"{gene_id}: produced {len(wt_frags_core)} fragments, exceeds max_fragments={max_frags_cfg}."
            )
        wt_frags_clone = backend.add_adapters(wt_frags_core, gg_cfg=gg_cfg)

        # Store WT fragment entries and seed pools with WT fragments
        if library_mode == "fragment_single_mutation":
            for f in wt_frags_clone:
                wt_id = f"{gene_id}|WT|frag{f.frag_index:02d}"
                wt_fragment_entries.append((wt_id, f.cloning_seq))

            # Add WT fragments to every pool (so each pool is self-contained for ordering)
            for pool_i in range(1, len(wt_frags_clone) + 1):
                pools.setdefault(pool_i, [])
                for f in wt_frags_clone:
                    wt_id = f"{gene_id}|WT|frag{f.frag_index:02d}"
                    pools[pool_i].append((wt_id, f.cloning_seq))

        # Pick generator based on mutation spec
        if mut_mode == "explicit":
            variants_iter = generate_explicit_variants(
                gene_id=gene_id,
                baseline_cds=baseline_cds,
                forbidden_motifs=forbidden,
                codon_table=table,
                explicit_mutations=mut_payload,  # type: ignore[arg-type]
                require_wt_match=True,
            )
        else:
            variants_iter = generate_saturation_variants(
                gene_id=gene_id,
                baseline_cds=baseline_cds,
                forbidden_motifs=forbidden,
                codon_table=table,
                include_wt=True,
                positions_to_mutate=mut_payload,  # type: ignore[arg-type]
            )

        for var in variants_iter:
            variant_id = f"{gene_id}|{var.mutation_label}"

            # Always write full-length variants for auditing
            variant_fasta_entries.append((variant_id, var.cds))
            variant_prot_entries.append((variant_id, var.protein))

            # Split variant + adapters
            var_frags_core = backend.split_sequence_by_cut_points(
                var.cds, cut_points=cut_points, enforce_frame=enforce_frame
            )
            var_frags_clone = backend.add_adapters(var_frags_core, gg_cfg=gg_cfg)

            # Map mutation to a fragment index by nt coordinate
            mut_nt = 3 * (var.pos1 - 1)
            mutated_frag_index: Optional[int] = None
            for core in var_frags_core:
                if core.start_nt <= mut_nt < core.end_nt:
                    mutated_frag_index = core.frag_index
                    break
            if mutated_frag_index is None:
                raise RuntimeError(f"{variant_id}: could not map mutation nt={mut_nt} to any fragment.")

            # Build core ORF for validation
            if library_mode == "fragment_single_mutation":
                # WT cores everywhere except mutated fragment core from var
                pieces = []
                for wt_core, var_core in zip(wt_frags_core, var_frags_core):
                    if wt_core.frag_index == mutated_frag_index:
                        pieces.append(var_core.core_seq)
                    else:
                        pieces.append(wt_core.core_seq)
                core_orf = "".join(pieces)

                # Save only the mutated fragment (cloning seq) for ordering
                mut_clone = next(f for f in var_frags_clone if f.frag_index == mutated_frag_index)
                mut_frag_id = f"{gene_id}|{var.mutation_label}|frag{mut_clone.frag_index:02d}"
                mut_fragment_entries.append((mut_frag_id, mut_clone.cloning_seq))

                # Add mutated fragment only to its pool
                pools[mutated_frag_index].append((mut_frag_id, mut_clone.cloning_seq))

                # Assembly plan row
                row = {
                    "gene_id": gene_id,
                    "construct_id": variant_id,
                    "mutation": var.mutation_label,
                    "mutated_frag_index": str(mutated_frag_index),
                }
                for f in wt_frags_clone:
                    if f.frag_index == mutated_frag_index:
                        row[f"frag{f.frag_index:02d}_id"] = mut_frag_id
                    else:
                        row[f"frag{f.frag_index:02d}_id"] = f"{gene_id}|WT|frag{f.frag_index:02d}"
                assembly_rows.append(row)

            else:
                # Classic behavior: validate ORF from the variant's own fragments
                core_orf = backend.reassemble_core_orf(var_frags_core)
                # Classic fragment output: all fragments per variant
                for f in var_frags_clone:
                    header = f"{variant_id}|frag{f.frag_index:02d}"
                    classic_fragment_entries.append((header, f.cloning_seq))

            # Validate assembled ORF
            vres = backend.check_construct(
                core_orf_seq=core_orf,
                codonopt_cfg=codonopt_cfg,
                final_forbidden_motifs=final_forbidden,
                require_no_internal_stops=bool(final_checks.get("require_no_internal_stops", True)),
                require_len_multiple_of_3=bool(final_checks.get("require_length_multiple_of_3", True)),
            )
            reassembled_entries.append((variant_id, core_orf))

            report_rows.append(
                {
                    "gene_id": gene_id,
                    "variant_id": variant_id,
                    "position_aa": str(var.pos1),
                    "wt_aa": var.wt_aa,
                    "target_aa": var.target_aa,
                    "mutation": var.mutation_label,
                    "codon_used": var.codon_used,
                    "mutated_frag_index": str(mutated_frag_index),
                    "passes_all": str(vres.passes_all),
                    "fail_reasons": ";".join(vres.fail_reasons),
                    "gc_fraction": str(vres.gc_fraction),
                    "max_homopolymer": str(vres.max_homopolymer),
                    "aa_len": str(vres.aa_len),
                    "n_frags": str(len(wt_frags_core)),
                    "library_mode": str(library_mode),
                    "mutation_mode": str(mut_mode),
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

    # --- WRITE OUTPUTS ---

    write_fasta(variants_fasta, variant_fasta_entries)
    write_fasta(variants_prot_fasta, variant_prot_entries)
    write_fasta(reassembled_fasta, reassembled_entries)

    # Classic fragments file only in classic mode
    if library_mode != "fragment_single_mutation":
        write_fasta(fragments_fasta, classic_fragment_entries)

    # TSV outputs
    def_hdr = [
        "gene_id",
        "variant_id",
        "position_aa",
        "wt_aa",
        "target_aa",
        "mutation",
        "codon_used",
        "mutated_frag_index",
        "passes_all",
        "fail_reasons",
        "gc_fraction",
        "max_homopolymer",
        "aa_len",
        "n_frags",
        "library_mode",
        "mutation_mode",
    ]
    write_tsv(report_tsv, report_rows, header=def_hdr)
    write_tsv(hits_tsv, hit_rows, header=["variant_id", "motif", "start", "end", "context"])

    # Pooling outputs
    if library_mode == "fragment_single_mutation":
        os.makedirs(pools_dir, exist_ok=True)

        write_fasta(wt_frags_fasta, wt_fragment_entries)
        write_fasta(mut_frags_fasta, mut_fragment_entries)

        # pool FASTAs
        for frag_index, entries in sorted(pools.items(), key=lambda kv: kv[0]):
            if frag_index > max_frags_cfg:
                continue
            pool_path = os.path.join(pools_dir, f"pool_frag{frag_index:02d}.fasta")
            # De-dup (header, seq) pairs while preserving order
            seen = set()
            uniq = []
            for hid, seq in entries:
                key = (hid, seq)
                if key not in seen:
                    uniq.append((hid, seq))
                    seen.add(key)
            write_fasta(pool_path, uniq)

        # assembly plan
        plan_header = ["gene_id", "construct_id", "mutation", "mutated_frag_index"] + [
            f"frag{i:02d}_id" for i in range(1, max_frags_cfg + 1)
        ]
        write_tsv(assembly_plan_tsv, assembly_rows, header=plan_header)
