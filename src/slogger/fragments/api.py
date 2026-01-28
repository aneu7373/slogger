from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Protocol, Sequence


@dataclass(frozen=True)
class FragmentCore:
    """Adapter-free fragment produced by splitting a CDS."""

    frag_index: int  # 1-based
    start: int       # 0-based
    end: int         # 0-based, exclusive
    core_seq: str


@dataclass(frozen=True)
class FragmentCloning(FragmentCore):
    """Fragment after Golden Gate adapters have been added."""

    cloning_seq: str


@dataclass(frozen=True)
class MotifHit:
    motif: str
    start: int  # 0-based
    end: int    # 0-based, exclusive
    context: str


@dataclass(frozen=True)
class ValidationResult:
    passes_all: bool
    fail_reasons: List[str]
    gc_fraction: float
    max_homopolymer: int
    aa_len: int
    hits: List[MotifHit]


class FragmentBackend(Protocol):
    """Stable SLOGGER-facing API for fragment/adapters/checks backends."""

    def split_sequence_by_cut_points(self, seq: str, cut_points: List[int], enforce_frame: bool = True) -> List[FragmentCore]:
        ...

    def add_adapters(self, fragments: List[FragmentCore], gg_cfg: List[Dict[str, str]]) -> List[FragmentCloning]:
        ...

    def reassemble_core_orf(self, fragments: List[FragmentCore]) -> str:
        ...

    def check_construct(
        self,
        core_orf_seq: str,
        codonopt_cfg: Dict[str, Any],
        final_forbidden_motifs: Optional[Sequence[str]] = None,
        require_no_internal_stops: bool = True,
        require_len_multiple_of_3: bool = True,
    ) -> ValidationResult:
        ...
