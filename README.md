# sirna-consensus

Consensus de novo siRNA design. Give it one transcript; get back ranked candidate
guides that **several independent open-source tools agreed on**.

```
FASTA  ->  fan out to N external tools  ->  normalise onto one axis
       ->  vote at the site level       ->  your filters  ->  rank  ->  report
```

This is not another design algorithm. It wraps tools that already exist, refuses
to let any single one of them decide, and keeps the survival rules in your hands.

---

## Why consensus, and why *site* level

Exact sequence intersection between three tools is usually **empty**: they emit
different core lengths (19 vs 21), different offsets (2 nt overhang conventions),
and different strands for the *same* duplex. So agreement is defined on the site:

> two candidates agree when they sit on the same transcript, start within
> `--offset-tolerance` (default 2 nt) of each other, and share all but 2 bases.

A site's vote is the number of **distinct tools** that proposed it — one verbose
tool cannot outvote the others.

Sequence filters are applied **after** voting. If the tools filtered first, "use
your own thresholds" would be a fiction.

## What this pipeline adds (and what it does not)

| It does | It does not |
|---|---|
| one coordinate axis for every tool's output | re-implement window enumeration |
| independent votes + explicit consensus tier | train yet another efficacy model |
| your filters, applied last | claim that agreement equals knockdown |
| full provenance (`run.json`, per-tool status) | hide a tool that failed to run |

The tools keep doing their own enumeration: they take the FASTA directly. The
only coordinate work here is *placing what comes back*.

---

## Install

```bash
python -m venv .venv && . .venv/bin/activate     # or: uv venv .venv
pip install -e ".[dev]"
sirna-consensus tools        # shows which external tools are usable
```

`sirna-consensus` itself has **no runtime dependencies**. The dependencies live
in each tool's own environment, so their pins never collide.

### Installing the external tools

Each tool runs in its own environment behind a process boundary. Put them
somewhere and point the environment at them.

**1. siRNAdesigner** — rule-based, fast, no external binaries.

```bash
git clone --depth 1 https://github.com/anazhmetdin/siRNAdesigner
python -m venv tools-venv && tools-venv/bin/pip install numpy biopython
export SIRNA_DESIGNER_SCRIPT=$PWD/siRNAdesigner/siRNADesigner.py
export SIRNA_CONSENSUS_PYTHON=$PWD/tools-venv/bin/python
```

**2. SIREN** — 21-mer off-target screen via RNAhybrid.

```bash
git clone --depth 1 https://github.com/pablovargasmejia/SIREN
tools-venv/bin/pip install ./SIREN/dist/siren_rnai-0.1.5-py3-none-any.whl
export SIREN_BIN=$PWD/tools-venv/bin/SIREN
# RNAhybrid is a C binary, not on PyPI. Without root:
mkdir -p toolbin && cd toolbin
apt-get download rnahybrid libg20
dpkg-deb -x rnahybrid_*.deb ext && dpkg-deb -x libg20_*.deb ext2
export SIRNA_CONSENSUS_TOOL_BIN=$PWD/ext/usr/bin
export SIRNA_CONSENSUS_TOOL_LIB=$PWD/ext2/usr/lib
```

**3. siRNAforge** — thermodynamics + off-target, from its pinned image.

```bash
docker pull ghcr.io/austin-s-h/sirnaforge:sha-4c05e4b
export SIRNAFORGE_IMAGE=ghcr.io/austin-s-h/sirnaforge:sha-4c05e4b
```

Run `sirna-consensus tools` until every tool you want reads `OK`.

### Environment variables

| variable | meaning |
|---|---|
| `SIRNA_CONSENSUS_PYTHON` | interpreter used for legacy tool scripts |
| `SIRNA_CONSENSUS_TOOL_BIN` | prepended to `PATH` for tool subprocesses (RNAhybrid) |
| `SIRNA_CONSENSUS_TOOL_LIB` | prepended to `LD_LIBRARY_PATH` for tool subprocesses |
| `SIRNA_CONSENSUS_TIMEOUT` | per-tool timeout in seconds (default 1800) |
| `SIRNA_CONSENSUS_TRANSCRIPTOME` | cDNA FASTA used for off-target search |
| `SIRNA_DESIGNER_SCRIPT` / `SIREN_BIN` / `SIRNAFORGE_IMAGE` | explicit tool locations |
| `SIRNAFORGE_SUBCMD` | `design` (default, fast) or `workflow` (adds off-target) |

---

## Use

```bash
sirna-consensus run \
  --fasta NM_001657.4.fa \
  --out out/ \
  --tools sirnadesigner,siren,sirnaforge \
  --min-votes 2 --top-k 10
```

Outputs in `--out`:

| file | contents |
|---|---|
| `candidates.tsv` | every consensus site, with votes, flags and per-tool scores |
| `top_guides.fa` | order-ready guide RNA, headers carry rank / tier / site / GC |
| `report.html` | standalone summary: tool status, normalisation, top table |
| `run.json` | machine-readable record: input hash, config, per-tool outcome |

Python API:

```python
from sirna_consensus import Config, run

result = run(Config(fasta="NM_001657.4.fa", out_dir="out/", min_votes=2))
for c in result["top"]:
    print(c.votes, c.core.transcript_id, c.core.start, c.core.guide)
```

### Reading the output

`candidates.tsv` is 1-based and inclusive on the input transcript. `*_dna`
columns align back to the transcript; `*_rna` columns are what you order. Guides
are the reverse complement of the target window — **no stage invents a base**.

`consensus_tier` is `unanimous` or `N-of-M`, where M is the number of tools that
actually ran. A tool that failed is recorded in `run.json` and reduces M rather
than silently vanishing.

---

## Filters

Hard (reject) — GC outside 30–55%, homopolymer ≥4, low sequence complexity,
guide 5′ end not A/U, internal palindrome.
Soft (annotate) — GU-rich guide, seed GC extremes, sense 3′ end.

These are the published design rules, deliberately rule-based so the ranking does
not collapse into "whichever model we happened to bundle".

## Limitations — read these before quoting a number

* **Consensus ≠ efficacy.** Three tools agreeing raises confidence in the site, it
  does not prove knockdown.
* **A tool that returns only its own Top-N is indistinguishable from one that
  filtered a site out.** Widen each tool's internal `top-n` where it allows it.
* **Tool scores are not comparable across tools** and are only used as a
  tie-break *inside* an equal-vote group.
* **Off-target scope is whatever the tool was given.** Pass
  `--transcriptome` to give them a shared reference, otherwise each tool screens
  against its own default (or nothing).
* **Consensus shrinks the list.** Unanimous sites are typically single digits — a
  feature if you want 5–10 orderable duplexes, a gap if you want exhaustive
  coverage. Singletons are kept by default and tagged by tier; drop them with
  `--exclude-singletons`.
* **Upstream bitrot is real.** siRNAdesigner needed a Biopython compatibility
  shim and an output flush; that shim is `tools/_shim.py` and leaves upstream
  untouched.

## Development

```bash
pip install -e ".[dev]"
pytest -q
```

The test suite is hermetic: it uses a fake tool, so CI needs no network, no
docker and no external binaries. `sirna-consensus tools` is the integration
check for a real machine.

## Licence

MIT for this package. Each wrapped tool keeps its own licence — check them
before redistribution (SIREN is GPL-3.0, siRNAforge is MIT).
