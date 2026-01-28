from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, List, Sequence, Tuple
from Bio.Seq import Seq

from .codon_table_xlsx import CodonTable, best_codons

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

@dataclass(frozen=True)
class Variant:
    gene_id: str
    pos1: int
    wt_aa: str
    target_aa: str
    mutation_label: str  # e.g. A123V or A123A
    cds: str
    protein: str
    codon_used: str  # codon at site (or "WT")

def _has_any_motif(seq: str, motifs: Sequence[str]) -> bool:
    s = seq.upper()
    for m in motifs:
        mm = str(m).upper().strip()
        if mm and mm in s:
            return True
    return False

def _swap_codon_try_synonyms(
    cds: str,
    pos1: int,
    target_aa: str,
    table: CodonTable,
    forbidden_motifs: Sequence[str],
) -> Tuple[str, str]:
    start = (pos1 - 1) * 3
    end = start + 3
    prefix = cds[:start]
    suffix = cds[end:]
    candidates = best_codons(table, target_aa)

    for codon in candidates:
        trial = prefix + codon + suffix
        if not _has_any_motif(trial, forbidden_motifs):
            return trial, codon

    # all fail: return top codon (validation will flag motifs later)
    return prefix + candidates[0] + suffix, candidates[0]

def generate_saturation_variants(
    gene_id: str,
    baseline_cds: str,
    forbidden_motifs: Sequence[str],
    codon_table: CodonTable,
    include_wt: bool = True,
) -> Iterable[Variant]:
    cds = baseline_cds.upper()
    prot = str(Seq(cds).translate(to_stop=False))
    if len(cds) != 3 * len(prot):
        raise ValueError(f"CDS/protein mismatch for {gene_id}")

    for pos1, wt in enumerate(prot, start=1):
        if wt == "*":
            continue
        for aa in AA_ORDER:
            if not include_wt and aa == wt:
                continue
            if aa == wt:
                mut_cds = cds
                used = "WT"
            else:
                mut_cds, used = _swap_codon_try_synonyms(
                    cds=cds, pos1=pos1, target_aa=aa, table=codon_table, forbidden_motifs=forbidden_motifs
                )
            mut_prot = prot[:pos1-1] + aa + prot[pos1:]
            label = f"{wt}{pos1}{aa}"
            yield Variant(
                gene_id=gene_id,
                pos1=pos1,
                wt_aa=wt,
                target_aa=aa,
                mutation_label=label,
                cds=mut_cds,
                protein=mut_prot,
                codon_used=used,
            )
