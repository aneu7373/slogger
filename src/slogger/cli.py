from __future__ import annotations
import argparse
import os
from .pipeline import run_pipeline

def main():
    p = argparse.ArgumentParser(prog="slogger", description="SLOGGER - saturation libraries, Golden-Gate excision ready.")
    sub = p.add_subparsers(dest="cmd", required=True)

    run = sub.add_parser("run", help="Run SLOGGER pipeline with a FROGGER-compatible YAML config.")
    run.add_argument("config", help="Path to YAML config.")
    run.add_argument("--outdir", default="out", help="Output directory.")

    args = p.parse_args()

    if args.cmd == "run":
        run_pipeline(args.config, args.outdir)
