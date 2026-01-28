from __future__ import annotations

from typing import Any, Dict, List, Optional, Sequence

from .api import FragmentCore, FragmentCloning, MotifHit, ValidationResult, FragmentBackend

from slogger.frogger_compat.v0_1 import split_sequence_by_cut_points as _split
from slogger.frogger_compat.v0_1 import add_adapters as _add
from slogger.frogger_compat.v0_1 import check_construct as _check


class FroggerBackendV0_1(FragmentBackend):
    """Backend that preserves FROGGER v0.1 fragment/adapters/check behavior."""

    def split_sequence_by_cut_points(self, seq: str, cut_points: List[int], enforce_frame: bool = True) -> List[FragmentCore]:
        frags = _split(seq=seq, cut_points=cut_points, enforce_frame=enforce_frame)
        out: List[FragmentCore] = []
        for idx, (a, b, core) in enumerate(frags, start=1):
            out.append(FragmentCore(frag_index=idx, start=a, end=b, core_seq=core))
        return out

    def add_adapters(self, fragments: List[FragmentCore], gg_cfg: List[Dict[str, str]]) -> List[FragmentCloning]:
        # Convert to FROGGER dict shape
        frogger_frags = [
            {
                "frag_index": f.frag_index,
                "start": f.start,
                "end": f.end,
                "core_seq": f.core_seq,
            }
            for f in fragments
        ]
        out_dicts = _add(frogger_frags, gg_cfg)
        # Convert back
        out: List[FragmentCloning] = []
        for d in out_dicts:
            out.append(
                FragmentCloning(
                    frag_index=int(d["frag_index"]),
                    start=int(d.get("start", 0)),
                    end=int(d.get("end", 0)),
                    core_seq=str(d["core_seq"]),
                    cloning_seq=str(d["cloning_seq"]),
                )
            )
        return out

    def reassemble_core_orf(self, fragments: List[FragmentCore]) -> str:
        # FROGGER semantics are simple concatenation in order
        frags_sorted = sorted(fragments, key=lambda f: f.frag_index)
        return "".join(f.core_seq for f in frags_sorted)

    def check_construct(
        self,
        core_orf_seq: str,
        codonopt_cfg: Dict[str, Any],
        final_forbidden_motifs: Optional[Sequence[str]] = None,
        require_no_internal_stops: bool = True,
        require_len_multiple_of_3: bool = True,
    ) -> ValidationResult:
        result = _check(
            construct={"core_orf_seq": core_orf_seq},
            codonopt_cfg=codonopt_cfg,
            final_forbidden_motifs=list(final_forbidden_motifs) if final_forbidden_motifs else None,
            require_no_internal_stops=require_no_internal_stops,
            require_len_multiple_of_3=require_len_multiple_of_3,
        )

        hits: List[MotifHit] = []
        for h in result.get("hits", []) or []:
            hits.append(
                MotifHit(
                    motif=str(h.get("motif", "")),
                    start=int(h.get("start", 0)),
                    end=int(h.get("end", 0)),
                    context=str(h.get("context", "")),
                )
            )

        return ValidationResult(
            passes_all=bool(result.get("passes_all")),
            fail_reasons=list(result.get("fail_reasons", []) or []),
            gc_fraction=float(result.get("gc_fraction", 0.0)),
            max_homopolymer=int(result.get("max_homopolymer", 0)),
            aa_len=int(result.get("aa_len", 0)),
            hits=hits,
        )
