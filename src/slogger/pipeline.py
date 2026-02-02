from __future__ import annotations

import glob
import os
import subprocess
from typing import Any, Dict, List, Optional, Tuple, Union

import yaml

from .io_fasta import read_fasta, write_fasta
from .codon_table_xlsx import load_codon_table_xlsx
from .mutate import (
    Variant,
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
    for c in candidates:
        if os.path.isdir(os.path.join(c, "codonopt")):
            existing = env.get("PYTHONPATH", "")
            env["PYTHONPATH"] = f"{c}{os.pathsep}{existing}" if existing else c
            break
    return env


def run_codonopt(cfg: Dict[str, Any], workdir: str) -> str:
    codonopt_cfg = cfg.get("codonopt", {}) or {}

    out_subdir = codonopt_cfg.get("out_subdir", "codonopt_out")
    out_dir = os.path.join(workdir, out_subdir)
    os.makedirs(out_dir, exist_ok=True)

    seq = codonopt_cfg.get("sequence")
    batch = codonopt_cfg.get("batch_table")
    if bool(seq) == bool(batch):
        raise ValueError("codonopt: provide exactly one of 'sequence' or 'batch_table'.")

    codon_table = codonopt_cfg.get("codon_table")
    if not codon_table:
        raise ValueError("codonopt.codon_table is required.")

    cmd: List[str] = ["python", "-m", "codonopt.main"]

    if seq:
        cmd += ["--sequence", str(seq)]
    if batch:
        cmd += ["--batch-table", str(batch)]

    cmd += ["--codon-table", str(codon_table)]

    if codonopt_cfg.get("codon_table_sheet"):
        cmd += ["--codon-table-sheet", str(codonopt_cfg["codon_table_sheet"])]

    cmd += ["--out", out_dir]

    if codonopt_cfg.get("verbose"):
        cmd += ["--verbose"]
    if codonopt_cfg.get("log_file"):
        cmd += ["--log-file", str(codonopt_cfg["log_file"])]

    for k, flag in [
        ("gc_min", "--gc-min"),
        ("gc_max", "--gc-max"),
        ("max_homopolymer", "--max-homopolymer"),
        ("min_codon_fraction", "--min-codon-fraction"),
        ("seed", "--seed"),
        ("n", "--n"),
    ]:
        if k in codonopt_cfg and codonopt_cfg[k] is not None:
            cmd += [flag, str(codonopt_cfg[k])]

    if codonopt_cfg.get("avoid_motifs"):
        cmd += ["--avoid-motifs", "|".join(map(str, codonopt_cfg["avoid_motifs"]))]

    env = _inject_codonopt_pythonpath(dict(os.environ), workdir)

    proc = subprocess.run(cmd, cwd=workdir, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            "codonopt failed:\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{proc.stdout}\n"
            f"STDERR:\n{proc.stderr}"
        )

    return _find_codonopt_output_fasta(out_dir)


MutationSpec = Tuple[str, Union[List[int], List[str]]]


def _parse_mutation_spec(cfg: Dict[str, Any], baseline_protein: str) -> MutationSpec:
    aa_len = len(baseline_protein)
    mut = cfg.get("mutation")

    if not mut:
        return ("saturation_positions", list(range(1, aa_len + 1)))

    if not isinstance(mut, dict):
        raise ValueError("mutation must be a mapping if provided.")

    if "explicit" in mut and mut["explicit"] is not None:
        return ("explicit", mut["explicit"])

    if "aa_range" in mut and mut["aa_range"] is not None:
        start, end = mut["aa_range"]
        return ("saturation_positions", list(range(int(start), int(end) + 1)))

    if "aa_positions" in mut and mut["aa_positions"] is not None:
        return ("saturation_positions", [int(x) for x in mut["aa_positions"]])

    return ("saturation_positions", list(range(1, aa_len + 1)))


def _normalize_protein(seq: str) -> str:
    s = (seq or "").strip().upper().replace(" ", "").replace("\n", "").replace("\r", "")
    if s.endswith("*"):
        s = s[:-1]
    return s


def _write_tsv(path: str, rows: List[Dict[str, Any]], header: List[str]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            out = []
            for h in header:
                v = r.get(h, "")
                out.append("" if v is None else str(v))
            f.write("\t".join(out) + "\n")


def _base_id(record_id: str) -> str:
    """
    codonopt may emit ids like:
      <original>|job0001|rep001
    Normalize back to <original> for joining with input FASTA IDs.
    """
    rid = (record_id or "").strip()
    if "|job" in rid:
        rid = rid.split("|job", 1)[0]
    # If someone has only |rep without |job, handle that too
    if "|rep" in rid:
        rid = rid.split("|rep", 1)[0]
    return rid


def run_pipeline(config_path: str, outdir: str) -> None:
    cfg = load_config(config_path)
    os.makedirs(outdir, exist_ok=True)

    workdir = os.path.dirname(os.path.abspath(config_path)) or "."

    inputs = cfg.get("inputs", {}) or {}
    input_fasta = inputs.get("fasta")
    input_type = str(inputs.get("input_type", "protein")).strip().lower()

    if input_type not in ("protein", "cds"):
        raise ValueError("inputs.input_type must be 'protein' or 'cds'.")

    if not input_fasta:
        raise ValueError("inputs.fasta is required.")

    input_fasta_path = os.path.join(workdir, input_fasta)

    codonopt_cfg = cfg.get("codonopt", {}) or {}
    codonopt_enable = bool(codonopt_cfg.get("enable", True))

    gg_cfg = (((cfg.get("golden_gate", {}) or {}).get("fragments", None)) or None)
    if gg_cfg is None:
        raise ValueError("golden_gate.fragments is required (adapters per fragment).")

    library_mode = ((cfg.get("library_mode", {}) or {}).get("type", "fragment_single_mutation") or "").strip()
    if library_mode != "fragment_single_mutation":
        raise ValueError("This build currently supports library_mode.type = fragment_single_mutation only.")

    # Read original input records for AA-preservation validation
    original_records = read_fasta(input_fasta_path)
    original_by_id = {r.id: r.seq for r in original_records}
    original_by_base = {_base_id(r.id): r.seq for r in original_records}

    # Baseline generation
    if codonopt_enable:
        if not codonopt_cfg.get("sequence") and not codonopt_cfg.get("batch_table"):
            cfg.setdefault("codonopt", {})
            cfg["codonopt"]["sequence"] = input_fasta

        optimized_fasta = run_codonopt(cfg, workdir)
        optimized = read_fasta(optimized_fasta)
    else:
        if input_type != "cds":
            raise ValueError("codonopt.enable=false is only allowed when inputs.input_type='cds'.")
        optimized = read_fasta(input_fasta_path)

    # Validate AA preservation for baseline optimized sequences
    for rec in optimized:
        base = _base_id(rec.id)
        if base not in original_by_base:
            raise ValueError(
                f"Record id '{rec.id}' (base='{base}') is present in codonopt output but not in input FASTA. "
                "IDs must match (or share a common base prefix before '|job'/'|rep') for AA-preservation validation."
            )

        optimized_prot = _normalize_protein(translate_cds(rec.seq))
        if input_type == "protein":
            expected = _normalize_protein(original_by_base[base])
        else:
            expected = _normalize_protein(translate_cds(original_by_base[base]))

        if optimized_prot != expected:
            raise ValueError(
                f"AA sequence mismatch after baseline optimization for record '{rec.id}' (base='{base}').\n"
                f"Expected (from input): {expected[:120]}{'...' if len(expected) > 120 else ''}\n"
                f"Got (from optimized CDS): {optimized_prot[:120]}{'...' if len(optimized_prot) > 120 else ''}\n"
                "This indicates either an input frame/translation issue or codonopt altered the amino-acid sequence."
            )

    # Codon table used for mutation codon selection
    if not codonopt_cfg.get("codon_table"):
        raise ValueError("codonopt.codon_table is required (used for mutation-site codon selection).")

    table = load_codon_table_xlsx(
        codonopt_cfg["codon_table"],
        codonopt_cfg.get("codon_table_sheet"),
    )

    forbidden = _collect_forbidden_motifs(cfg)
    backend = FroggerBackendV0_1()

    # Outputs
    baseline_out = os.path.join(outdir, "baseline_optimized.fasta")
    variants_out = os.path.join(outdir, "variants.fasta")

    write_fasta(baseline_out, [(r.id, r.seq) for r in optimized])

    out_cfg = cfg.get("outputs", {}) or {}
    wt_frags_path = os.path.join(outdir, out_cfg.get("wt_fragments_fasta", "wt_fragments_with_gg.fasta"))
    mut_frags_path = os.path.join(outdir, out_cfg.get("mut_fragments_fasta", "mut_fragments_with_gg.fasta"))
    report_path = os.path.join(outdir, out_cfg.get("report_tsv", "report.tsv"))
    assembly_plan_path = os.path.join(outdir, out_cfg.get("assembly_plan_tsv", "assembly_plan.tsv"))
    pools_dir = os.path.join(outdir, out_cfg.get("pools_dir", "pools"))
    os.makedirs(pools_dir, exist_ok=True)

    # Collections
    variant_entries: List[Tuple[str, str]] = []
    wt_fragment_entries: List[Tuple[str, str]] = []
    mut_fragment_entries: List[Tuple[str, str]] = []

    pools: Dict[int, List[Tuple[str, str]]] = {}
    report_rows: List[Dict[str, Any]] = []
    assembly_rows: List[Dict[str, Any]] = []

    for rec in optimized:
        gene_id_full = rec.id
        gene_id = _base_id(gene_id_full)

        baseline_cds = rec.seq
        baseline_prot = translate_cds(baseline_cds)

        mut_mode, mut_payload = _parse_mutation_spec(cfg, baseline_prot)
        cut_points, enforce_frame = _merge_cut_points(cfg, gene_id)

        wt_frags_core = backend.split_sequence_by_cut_points(
            baseline_cds, cut_points=cut_points, enforce_frame=enforce_frame
        )
        wt_frags_clone = backend.add_adapters(wt_frags_core, gg_cfg=gg_cfg)

        n_frags = len(wt_frags_clone)
        if n_frags < 1:
            raise RuntimeError(f"{gene_id}: no fragments produced. Check cut points.")

        # WT fragments + seed pools with WT
        for f in wt_frags_clone:
            wt_id = f"{gene_id}|WT|frag{f.frag_index:02d}"
            wt_fragment_entries.append((wt_id, f.cloning_seq))

        for i in range(1, n_frags + 1):
            pools.setdefault(i, [])
            for f in wt_frags_clone:
                pools[i].append((f"{gene_id}|WT|frag{f.frag_index:02d}", f.cloning_seq))

        # Variants
        if mut_mode == "explicit":
            variants_iter = generate_explicit_variants(
                gene_id=gene_id,
                baseline_cds=baseline_cds,
                forbidden_motifs=forbidden,
                codon_table=table,
                explicit_mutations=mut_payload,  # type: ignore[arg-type]
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
            assert isinstance(var, Variant)

            variant_id = f"{gene_id}|{var.mutation_label}"
            variant_entries.append((variant_id, var.cds))

            # Map mutation position -> fragment index by cumulative lengths
            mut_nt = 3 * (var.pos1 - 1)
            offset = 0
            mutated_frag_index: Optional[int] = None
            for core in wt_frags_core:
                length = len(core.core_seq)
                if offset <= mut_nt < offset + length:
                    mutated_frag_index = core.frag_index
                    break
                offset += length

            if mutated_frag_index is None:
                raise RuntimeError(
                    f"{gene_id}|{var.mutation_label}: failed to map mutation nt={mut_nt} to a fragment "
                    f"(total_len={offset}, n_frags={len(wt_frags_core)})."
                )

            # Mutant fragments for this variant
            var_frags_core = backend.split_sequence_by_cut_points(
                var.cds, cut_points=cut_points, enforce_frame=enforce_frame
            )
            var_frags_clone = backend.add_adapters(var_frags_core, gg_cfg=gg_cfg)

            mut_frag = next(f for f in var_frags_clone if f.frag_index == mutated_frag_index)
            mut_frag_id = f"{gene_id}|{var.mutation_label}|frag{mutated_frag_index:02d}"

            mut_fragment_entries.append((mut_frag_id, mut_frag.cloning_seq))
            pools[mutated_frag_index].append((mut_frag_id, mut_frag.cloning_seq))

            report_rows.append(
                {
                    "gene_id": gene_id,
                    "variant_id": variant_id,
                    "mutation": var.mutation_label,
                    "position_aa": var.pos1,
                    "wt_aa": var.wt_aa,
                    "target_aa": var.target_aa,
                    "codon_used": var.codon_used,
                    "mutated_frag_index": mutated_frag_index,
                    "n_frags": n_frags,
                }
            )

            row: Dict[str, Any] = {
                "gene_id": gene_id,
                "variant_id": variant_id,
                "mutation": var.mutation_label,
                "mutated_frag_index": mutated_frag_index,
            }
            for f in wt_frags_clone:
                frag_key = f"frag{f.frag_index:02d}_id"
                if f.frag_index == mutated_frag_index:
                    row[frag_key] = mut_frag_id
                else:
                    row[frag_key] = f"{gene_id}|WT|frag{f.frag_index:02d}"
            assembly_rows.append(row)

    # Write FASTA outputs
    write_fasta(variants_out, variant_entries)
    write_fasta(wt_frags_path, wt_fragment_entries)
    write_fasta(mut_frags_path, mut_fragment_entries)

    # Write pools
    for idx, entries in sorted(pools.items(), key=lambda kv: kv[0]):
        seen = set()
        uniq = []
        for h, s in entries:
            key = (h, s)
            if key not in seen:
                uniq.append((h, s))
                seen.add(key)
        write_fasta(os.path.join(pools_dir, f"pool_frag{idx:02d}.fasta"), uniq)

    # Write report.tsv
    report_header = [
        "gene_id",
        "variant_id",
        "mutation",
        "position_aa",
        "wt_aa",
        "target_aa",
        "codon_used",
        "mutated_frag_index",
        "n_frags",
    ]
    _write_tsv(report_path, report_rows, report_header)

    # Write assembly_plan.tsv
    max_idx = 0
    for r in assembly_rows:
        for k in r.keys():
            if k.startswith("frag") and k.endswith("_id"):
                try:
                    num = int(k[4:6])
                    max_idx = max(max_idx, num)
                except Exception:
                    pass

    plan_header = ["gene_id", "variant_id", "mutation", "mutated_frag_index"] + [
        f"frag{i:02d}_id" for i in range(1, max_idx + 1)
    ]
    _write_tsv(assembly_plan_path, assembly_rows, plan_header)
