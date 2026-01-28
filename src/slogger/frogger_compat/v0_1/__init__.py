"""FROGGER compatibility layer (ported snapshot).

Contains ported, mostly-verbatim logic from FROGGER:
- splitters.py
- gg_adapters.py
- checks.py

SLOGGER should call these through slogger.fragments.frogger_backend, not directly.
"""

from .splitters import split_sequence_by_cut_points
from .gg_adapters import add_adapters
from .checks import check_construct

__all__ = [
    "split_sequence_by_cut_points",
    "add_adapters",
    "check_construct",
]
