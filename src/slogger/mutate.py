from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_SET = set(AA_ORDER)

# Standard genetic code (DNA codons)
CODON_TO_AA = {
    # Phenylalanine / Leucine
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Isoleucine / Methionine
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    # Valine
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Serine / Proline / Threonine / Alanine
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine / Histidine / Glutamine / Asparagine / Lysine / Aspartate / Glutamate
    "TAT": "Y", "TAC": "Y",
    "CAT": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    # Cysteine / Tryptophan / Arginine / Glycine
    "TGT": "C", "TGC": "C", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stops
    "TAA": "*", "TAG": "*", "TGA": "*",
}

MUT_RE = re.compile(r"^([A-Za-z])(\d+)([A-Za-z])$")


@dataclass(frozen=True)
class Variant:
    gene_id: str
    pos1: int
    wt_aa: str
    target_aa: str
    mutation_label: str  # e.g. A123V or A123A
    codon_used: str      # codon inserted at mutated site, or "WT"
    protein: str
    cds: str


def _normalize_dna(seq: str) -> str:
    return seq.upper().replace("U", "T").replace(" ", "").replace("\n", "").replace("\r", "")


def translate_cds(cds: str) -> str:
    cds = _normalize_dna(cds)
    if len(cds) % 3 != 0:
        raise ValueError("CDS length is not a multiple of 3.")
    prot = []
    for i in range(0, len(cds), 3):
        codon = cds[i:i + 3]
        aa = CODON_TO_AA.get(codon)
        if aa is None:
            raise ValueError(f"Invalid/unknown codon: {codon}")
        prot.append(aa)
    return "".join(prot)


def _contains_any_motif(seq: str, motifs: Sequence[str]) -> bool:
    s = seq.upper()
    for m in motifs:
        mm = str(m).upper()
        if mm and mm in s:
            return True
    return False


def _codon_choices_for_aa(
    target_aa: str,
    codon_table: Dict[str, List[Tuple[str, float]]],
) -> List[str]:
    """
    codon_table expected shape:
      {"A": [("GCT", 0.35), ("GCC", 0.28), ...], "R": [...], ...}
    Sorted high->low frequency preferred.
    """
    aa = target_aa.upper()
    if aa not in codon_table or not codon_table[aa]:
        raise ValueError(f"No codon choices provided for amino acid '{aa}' in codon_table.")
    out: List[str] = []
    for cod, _freq in codon_table[aa]:
        out.append(_normalize_dna(cod))
    # de-dup while preserving order
    seen = set()
    uniq = []
    for c in out:
        if c not in seen:
            uniq.append(c)
            seen.add(c)
    return uniq


def choose_mutation_codon(
    baseline_cds: str,
    codon_index0: int,
    target_aa: str,
    forbidden_motifs: Sequence[str],
    codon_table: Dict[str, List[Tuple[str, float]]],
) -> Tuple[str, str]:
    """
    Choose the highest-frequency codon for target_aa unless it introduces a forbidden motif.
    Deterministically tries synonymous codons in descending frequency order.

    Returns: (codon_used, mutated_cds)
    """
    baseline_cds = _normalize_dna(baseline_cds)
    start = 3 * codon_index0
    end = start + 3
    if start < 0 or end > len(baseline_cds):
        raise IndexError("codon_index0 is out of bounds for baseline_cds.")

    choices = _codon_choices_for_aa(target_aa, codon_table)
    for codon in choices:
        candidate = baseline_cds[:start] + codon + baseline_cds[end:]
        if not _contains_any_motif(candidate, forbidden_motifs):
            return codon, candidate

    # If none work, return the top choice but let downstream validation flag it.
    top = choices[0]
    candidate = baseline_cds[:start] + top + baseline_cds[end:]
    return top, candidate


def parse_explicit_mutation_string(mut: str) -> Tuple[str, int, str]:
    """
    Parse a mutation string like 'T24V' into (wt_aa, pos1, target_aa).
    Allows lowercase input. Validates AAs are in the standard 20.
    """
    m = MUT_RE.match(mut.strip())
    if not m:
        raise ValueError(f"Invalid mutation string '{mut}'. Expected format like 'T24V'.")

    wt = m.group(1).upper()
    pos1 = int(m.group(2))
    tgt = m.group(3).upper()

    if wt not in AA_SET:
        raise ValueError(f"Invalid WT amino acid '{wt}' in '{mut}'. Must be one of {AA_ORDER}.")
    if tgt not in AA_SET:
        raise ValueError(f"Invalid target amino acid '{tgt}' in '{mut}'. Must be one of {AA_ORDER}.")
    if pos1 < 1:
        raise ValueError(f"Invalid position {pos1} in '{mut}'. Positions are 1-based and must be >= 1.")

    return wt, pos1, tgt


def generate_saturation_variants(
    *,
    gene_id: str,
    baseline_cds: str,
    forbidden_motifs: Sequence[str],
    codon_table: Dict[str, List[Tuple[str, float]]],
    include_wt: bool = True,
    positions_to_mutate: Optional[Sequence[int]] = None,  # 1-based AA indices
) -> Iterable[Variant]:
    """
    Generate exhaustive single-site saturation variants across requested AA positions.
    Default: all positions.

    Includes WT per position (e.g., A123A) if include_wt=True.
    """
    baseline_cds = _normalize_dna(baseline_cds)
    baseline_prot = translate_cds(baseline_cds)

    aa_len = len(baseline_prot)
    if positions_to_mutate is None:
        positions = list(range(1, aa_len + 1))
    else:
        positions = list(positions_to_mutate)

    # Validate positions
    for p in positions:
        if p < 1 or p > aa_len:
            raise ValueError(f"Requested mutation position {p} outside protein length {aa_len} for {gene_id}.")

    for pos1 in positions:
        wt = baseline_prot[pos1 - 1]
        codon_index0 = pos1 - 1

        if include_wt:
            mut_label = f"{wt}{pos1}{wt}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=wt,
                target_aa=wt,
                mutation_label=mut_label,
                codon_used="WT",
                protein=baseline_prot,
                cds=baseline_cds,
            )

        for target in AA_ORDER:
            if target == wt:
                continue
            codon_used, mutated_cds = choose_mutation_codon(
                baseline_cds=baseline_cds,
                codon_index0=codon_index0,
                target_aa=target,
                forbidden_motifs=forbidden_motifs,
                codon_table=codon_table,
            )
            mutated_prot = translate_cds(mutated_cds)
            mut_label = f"{wt}{pos1}{target}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=wt,
                target_aa=target,
                mutation_label=mut_label,
                codon_used=codon_used,
                protein=mutated_prot,
                cds=mutated_cds,
            )


def generate_explicit_variants(
    *,
    gene_id: str,
    baseline_cds: str,
    forbidden_motifs: Sequence[str],
    codon_table: Dict[str, List[Tuple[str, float]]],
    explicit_mutations: Sequence[str],  # strings like "T24V"
    require_wt_match: bool = True,
) -> Iterable[Variant]:
    """
    Generate exactly the variants specified in explicit_mutations.

    - If require_wt_match=True, we verify the WT letter matches the baseline protein at that position.
      This prevents off-by-one and wrong-sequence mistakes.
    - WT entries are allowed (e.g. T24T) and will yield the baseline CDS unchanged.
    """
    baseline_cds = _normalize_dna(baseline_cds)
    baseline_prot = translate_cds(baseline_cds)
    aa_len = len(baseline_prot)

    # de-dup while preserving order
    seen = set()
    muts: List[str] = []
    for m in explicit_mutations:
        mm = str(m).strip()
        if mm and mm not in seen:
            muts.append(mm)
            seen.add(mm)

    for mstr in muts:
        wt, pos1, tgt = parse_explicit_mutation_string(mstr)
        if pos1 > aa_len:
            raise ValueError(f"{gene_id}: explicit mutation {mstr} exceeds protein length {aa_len}.")

        baseline_wt = baseline_prot[pos1 - 1]
        if require_wt_match and baseline_wt != wt:
            raise ValueError(
                f"{gene_id}: explicit mutation {mstr} WT='{wt}' does not match baseline AA at position {pos1} "
                f"(baseline has '{baseline_wt}')."
            )

        codon_index0 = pos1 - 1

        # WT "mutation"
        if tgt == baseline_wt:
            mut_label = f"{baseline_wt}{pos1}{baseline_wt}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=baseline_wt,
                target_aa=baseline_wt,
                mutation_label=mut_label,
                codon_used="WT",
                protein=baseline_prot,
                cds=baseline_cds,
            )
            continue

        codon_used, mutated_cds = choose_mutation_codon(
            baseline_cds=baseline_cds,
            codon_index0=codon_index0,
            target_aa=tgt,
            forbidden_motifs=forbidden_motifs,
            codon_table=codon_table,
        )
        mutated_prot = translate_cds(mutated_cds)
        mut_label = f"{baseline_wt}{pos1}{tgt}"
        yield Variant(
            gene_id=gene_id,
            pos1=pos1,
            wt_aa=baseline_wt,
            target_aa=tgt,
            mutation_label=mut_label,
            codon_used=codon_used,
            protein=mutated_prot,
            cds=mutated_cds,
        )
