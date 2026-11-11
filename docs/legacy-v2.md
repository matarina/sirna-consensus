# Legacy schema-v2 compatibility

This document records compatibility semantics for older `run.json` and `select`
outputs. It describes parsing and molecule identity only; the repository no longer
ships the historical ranker, training data, or benchmark artifacts.

## Pipeline behavior

`run` executes registered external tools and stores every normalized hit. `select`
replays those exact molecules or an explicit candidate TSV. The package never
enumerates a transcript locally. Generator tools provide candidates; SIREN remains an
evaluator and cannot create a design vote.

The default generators are siRNADesigner and siRNAforge. Their source/image pins,
versions, native scores, chemistry, coordinates, strands, and raw-output hashes are
recorded independently. Invalid coordinates, inconsistent strands, unknown chemistry,
and unsupported modifications are rejected.

## Molecule identity

Coordinates are 1-based inclusive on the supplied transcript. Binding cores,
complete guide/passenger strands, and documented overhangs are separate fields.
siRNAdesigner's staggered 19-bp duplex and siRNAforge's blunt 21-mer remain distinct
formats. A candidate is never truncated to fit another format.

Molecular IDs include strand and overhang identity. Occurrence keys include transcript
and position. Passenger-known and passenger-unknown records remain separate. Repeated
sites at different coordinates do not share evidence.

## Selection

Native ranking uses the scores emitted by each generator, with deterministic tie
breaking. Optional `model`, `fusion`, `support`, and `quality` strategies are available
only when the caller supplies an external scorer manifest. The adapter performs no
training and validates manifest, checkpoint, command, and artifact hashes before use.

```bash
sirna-consensus select --fasta target.fa --run native-run/run.json \
  --scorer-manifest /path/to/scorer_manifest.json \
  --strategy model --out scored
```

A scorer checkpoint and its training provenance are external inputs and are recorded in
`selection.json`; they are not bundled with this repository. Without a scorer manifest,
use `--strategy native` or an explicit candidate TSV.

## Specificity screening

A pinned Bowtie reference manifest can be supplied to `select`:

```bash
sirna-consensus select --fasta target.fa --run native-run/run.json --out screened \
  --screen-manifest /path/to/screen_manifest.json \
  --target-gene ENSG00000100644 --bowtie /path/to/bowtie
```

The screen checks paired binding regions on both strands. Exact or one-mismatch hits
in unintended genes exclude a guide from the screened panel; two-mismatch and seed
burdens are advisory. Missing mappings, failed execution, hash mismatches, missing
intended hits, and incomplete UTR projection leave a candidate unverified.

## Historical files

Schema-v1 records and old reports may be inspected from external archives, but they are
not part of this package. Reproducible efficacy claims require exact measured duplexes
and an independently maintained assay dataset. An unmeasured generated molecule is
unknown, never inactive.

Normal tests are hermetic. Adapter tests that invoke external tools require
`SIRNA_CONSENSUS_EXTERNAL_TESTS=1`.
