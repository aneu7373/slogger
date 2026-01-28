from typing import List, Tuple


def split_sequence_by_cut_points(
    seq: str,
    cut_points: List[int],
    enforce_frame: bool = True,
) -> List[Tuple[int, int, str]]:
    """
    cut_points: 0-based cut points into seq
    Example: [300,600] => [0:300],[300:600],[600:end]
    """
    if any(cp < 0 or cp > len(seq) for cp in cut_points):
        raise ValueError("Cut points must be within [0, len(seq)].")
    if sorted(cut_points) != list(cut_points):
        raise ValueError("Cut points must be sorted ascending.")
    if enforce_frame and any((cp % 3) != 0 for cp in cut_points):
        raise ValueError("Cut points must be multiples of 3 when enforce_frame=true.")

    points = [0] + list(cut_points) + [len(seq)]
    frags = []
    for a, b in zip(points[:-1], points[1:]):
        if b <= a:
            raise ValueError("Invalid cut points (non-increasing).")
        frags.append((a, b, seq[a:b]))
    return frags
