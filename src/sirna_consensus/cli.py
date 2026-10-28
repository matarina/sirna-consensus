"""Command line interface."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__, filters
from .pipeline import Config, run
from .tools import DEFAULT_TOOLS, available_tools, survey


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sirna-consensus",
        description="Consensus de novo siRNA design from several external tools.")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)

    run_parser = sub.add_parser("run", help="design candidates for one transcript")
    run_parser.add_argument("--fasta", required=True,
                            help="single-record transcript FASTA (the coordinate axis)")
    run_parser.add_argument("--out", required=True, help="output directory")
    run_parser.add_argument("--tools", default=",".join(DEFAULT_TOOLS),
                            help=f"comma separated; known: {','.join(available_tools())}")
    run_parser.add_argument("--min-votes", type=int, default=2,
                            help="generator tools that must agree (default 2); "
                                 "also enforced for the exported Top-K")
    run_parser.add_argument("--top-k", type=int, default=10)
    run_parser.add_argument("--min-spacing", type=int, default=30,
                            help="minimum nt between selected top candidates")
    run_parser.add_argument("--lengths", default="19,20,21")
    run_parser.add_argument("--offset-tolerance", type=int, default=2,
                            help="max nt offset for cores to be grouped into a site")
    run_parser.add_argument("--site-mode", choices=["bounded", "exact", "fuzzy"],
                            default="bounded",
                            help="bounded = anchored offset clustering (default); "
                                 "exact = identical start+length; fuzzy = single-linkage "
                                 "overlap, which can chain a dense transcript")
    run_parser.add_argument("--transcriptome", default="",
                            help="cDNA FASTA for tools that support off-target search")
    run_parser.add_argument("--jobs", type=int, default=3)
    run_parser.add_argument("--exclude-singletons", action="store_true",
                            help="drop singleton sites from candidates.tsv; "
                                 "Top-K already enforces --min-votes")
    run_parser.add_argument("--strict-specificity", action="store_true",
                            help="exclude candidates lacking a clean, "
                                 "transcriptome-backed specificity screen")
    run_parser.add_argument("--siren-max-offtargets", type=int, default=0,
                            help="max off-targets for a SIREN site to be clean evidence")
    run_parser.add_argument("--gc-min", type=float, default=filters.GC_MIN)
    run_parser.add_argument("--gc-max", type=float, default=filters.GC_MAX)
    run_parser.add_argument("--max-homopolymer", type=int, default=4)
    run_parser.add_argument("--gu-rich-cutoff", type=float, default=filters.GU_RICH_CUTOFF)
    run_parser.add_argument("--seed-gc-high", type=float, default=0.75)
    run_parser.add_argument("--seed-gc-low", type=float, default=0.25)

    sub.add_parser("tools", help="list registered tools and whether they are usable")
    return parser


def _cmd_tools() -> int:
    rows = survey()
    width = max(len(name) for name, _, _ in rows)
    for name, usable, detail in rows:
        mark = "OK " if usable else "-- "
        print(f"{mark}{name.ljust(width)}  {detail}")
    if not any(usable for _, usable, _ in rows):
        print("\nno tool is usable; see README 'Installing the external tools'",
              file=sys.stderr)
        return 1
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    lengths = tuple(int(x) for x in args.lengths.split(",") if x.strip())
    if not lengths or any(length <= 0 for length in lengths):
        sys.stderr.write("--lengths must be a comma-separated list of positive integers\n")
        return 2

    if args.min_votes < 1:
        sys.stderr.write("--min-votes must be >= 1\n")
        return 2
    if args.top_k < 1:
        sys.stderr.write("--top-k must be >= 1\n")
        return 2
    if args.min_spacing < 0:
        sys.stderr.write("--min-spacing must be >= 0\n")
        return 2
    if args.offset_tolerance < 0:
        sys.stderr.write("--offset-tolerance must be >= 0\n")
        return 2
    if args.jobs < 1:
        sys.stderr.write("--jobs must be >= 1\n")
        return 2
    if args.siren_max_offtargets < 0:
        sys.stderr.write("--siren-max-offtargets must be >= 0\n")
        return 2

    try:
        filter_config = filters.FilterConfig(
            gc_min=args.gc_min,
            gc_max=args.gc_max,
            max_homopolymer=args.max_homopolymer,
            gu_rich_cutoff=args.gu_rich_cutoff,
            seed_gc_high=args.seed_gc_high,
            seed_gc_low=args.seed_gc_low,
        )
    except ValueError as exc:
        sys.stderr.write(f"invalid filter configuration: {exc}\n")
        return 2

    tool_names = tuple(t.strip() for t in args.tools.split(",") if t.strip())
    if not tool_names:
        sys.stderr.write("--tools must list at least one tool\n")
        return 2
    if args.transcriptome and "siren" in tool_names:
        sys.stderr.write(
            "note: a transcriptome-wide SIREN screen is expensive; "
            "SIREN_THREADS controls RNAhybrid concurrency\n")

    config = Config(
        fasta=args.fasta,
        out_dir=args.out,
        tools=tool_names,
        lengths=lengths,
        min_votes=args.min_votes,
        top_k=args.top_k,
        min_spacing=args.min_spacing,
        offset_tolerance=args.offset_tolerance,
        transcriptome=args.transcriptome,
        jobs=args.jobs,
        site_mode=args.site_mode,
        include_singletons=not args.exclude_singletons,
        strict_specificity=args.strict_specificity,
        siren_max_offtargets=args.siren_max_offtargets,
        filter_config=filter_config,
    )
    result = run(config)
    counts = result["counts"]
    print(f"transcript {result['transcript_id']} ({result['transcript_length']} nt)")
    for tool in result["tools"]:
        print(f"  {tool['tool']:<16} {tool['role']:<10} {tool['status']:<12} "
              f"{tool['hits']:>6} hits  {tool['detail']}")
    print(f"  placed {counts['placed']} · sites {counts['sites']} · "
          f"generator votes {counts['generator_hits']} · "
          f"consensus candidates {counts['consensus_candidates']} · "
          f"exploratory singletons {counts['exploratory_singletons']} · "
          f"unanimous {counts['consensus_unanimous']} · "
          f"below vote threshold {counts['below_min_votes']} · "
          f"rejected {counts['rejected_by_filters']} · "
          f"passing {counts['passing']} · top {counts['top']} · "
          f"clean specificity {counts['evaluator_clean']}")
    print(f"wrote {args.out}/candidates.tsv, top_guides.fa, all_candidates.fa, "
          f"report.html, run.json")
    print(json.dumps({"tool": "sirna-consensus", "counts": counts}))
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command == "tools":
        return _cmd_tools()
    if args.command == "run":
        return _cmd_run(args)
    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
