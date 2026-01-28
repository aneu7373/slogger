from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable, List, Tuple
from Bio import SeqIO

@dataclass
class FastaRecord:
    id: str
    description: str
    seq: str

def read_fasta(path: str) -> List[FastaRecord]:
    recs: List[FastaRecord] = []
    for r in SeqIO.parse(path, "fasta"):
        recs.append(FastaRecord(id=r.id, description=r.description, seq=str(r.seq)))
    if not recs:
        raise ValueError(f"No FASTA records found: {path}")
    return recs

def write_fasta(path: str, entries: Iterable[Tuple[str, str]]) -> None:
    with open(path, "w", encoding="utf-8") as f:
        for h, s in entries:
            f.write(f">{h}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")
