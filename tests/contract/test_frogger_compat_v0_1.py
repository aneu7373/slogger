import pytest

from slogger.fragments import FroggerBackendV0_1


def test_split_and_adapters_roundtrip():
    backend = FroggerBackendV0_1()
    cds = "ATGAAACCCGGG"  # 12 nt
    cut_points = [6]  # 0:6, 6:12 (frame-aligned)
    frags = backend.split_sequence_by_cut_points(cds, cut_points, enforce_frame=True)
    assert [f.core_seq for f in frags] == ["ATGAAA", "CCCGGG"]
    assert [f.frag_index for f in frags] == [1, 2]

    gg_cfg = [
        {"left_adapter": "L1", "right_adapter": "R1"},
        {"left_adapter": "L2", "right_adapter": "R2"},
    ]
    cloning = backend.add_adapters(frags, gg_cfg)
    assert cloning[0].cloning_seq == "L1ATGAAAR1"
    assert cloning[1].cloning_seq == "L2CCCGGGR2"

    core_orf = backend.reassemble_core_orf(frags)
    assert core_orf == cds


def test_frame_enforcement_raises():
    backend = FroggerBackendV0_1()
    cds = "ATGAAACCCGGG"
    with pytest.raises(ValueError):
        backend.split_sequence_by_cut_points(cds, [5], enforce_frame=True)


def test_check_construct_motif_and_metrics():
    backend = FroggerBackendV0_1()

    # Contains AAAA motif twice overlapping positions: ATGAAAAAA...
    core_orf = "ATGAAAAAATGA"  # 12 nt, translates to MKK*
    codonopt_cfg = {
        "max_homopolymer": 5,
        "gc_min": 0.0,
        "gc_max": 1.0,
        "avoid_motifs": ["AAAA"],
    }
    res = backend.check_construct(
        core_orf_seq=core_orf,
        codonopt_cfg=codonopt_cfg,
        final_forbidden_motifs=[],
        require_no_internal_stops=True,
        require_len_multiple_of_3=True,
    )
    assert res.passes_all is False
    assert "forbidden_motif_or_site" in res.fail_reasons
    assert len(res.hits) >= 1
    assert res.max_homopolymer >= 4
