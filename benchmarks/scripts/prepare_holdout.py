"""Build the weight-independent holdout dataset from the audited snapshot.

The cohort (studies sharing no molecule with checkpoint training or the development
cohort) is read from ``cohort_manifest.json``; nothing is re-selected here.  Each
record's reported antisense 21-mer and 19-mer sense are placed on the publisher
accession sequence by exact search.  Coordinates are the publisher's cloned/reported
constructs, never current full transcripts.  Records that cannot be placed exactly
are excluded with a reason, not guessed.

Endpoint is ``percent_inhibition`` (threshold 70); the snapshot's Huesken-style
per-cent normalisation is a separate unit and is not used here.
"""
import argparse
import collections
import csv
import json
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from audit_exposure import parse_snapshot  # noqa: E402

from sirna_consensus.fasta import sha256, sha256_file  # noqa: E402
from sirna_consensus.types import Core, reverse_complement, to_dna  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def read_fasta(path: Path) -> dict[str, str]:
    """Reference sequences, with ambiguous bases replaced in place by 'A'.

    The library's reader requires a pure ACGT alphabet; the publisher's Lamin A
    construct (only) carries N runs. Replacement preserves 1-based coordinates, so
    any guide whose sense window is N-free still places exactly.
    """
    sequences, current, parts = {}, None, []
    for line in [*Path(path).read_text().splitlines(), ">END"]:
        if line.startswith(">"):
            if current:
                sequences[current] = re.sub(r"[^ACGT]", "A", to_dna("".join(parts)))
            current, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    return sequences


def safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def build(snapshot: Path, fasta: Path, cohort_path: Path) -> tuple[list, list, dict, list, dict]:
    if not cohort_path.is_file():
        raise FileNotFoundError(f'cohort manifest not found: {cohort_path}; run audit_exposure.py')
    holdout_studies = set(json.loads(cohort_path.read_text())["holdout_studies"])
    sequences = read_fasta(fasta)
    by_base = {}
    for key in sequences:
        by_base.setdefault(key.split(".")[0], key)
    records, candidates, exclusions, transcripts = [], [], [], {}
    by_key = defaultdict(list)
    def exclude(record_id, study, reason):
        exclusions.append({"record_id": record_id, "study": study, "reason": reason})

    for row in parse_snapshot(snapshot):
        study = row["Authors"]
        if study not in holdout_studies:
            continue
        rid = f"snapshot:{row['ID']}"
        reported = (row.get("Accession_number", "") or "").strip()
        base = reported.split(".")[0]
        accession = base if base in sequences else by_base.get(base)
        if not accession:
            exclude(rid, study, f"accession {reported!r} not in reference")
            continue
        tx = sequences[accession]
        sense = to_dna(row.get("Sense_19mer", ""))
        guide = to_dna(row.get("Antisense_21mer", ""))
        if len(sense) != 19 or len(guide) not in (19, 21) or reverse_complement(guide[:19]) != sense:
            exclude(rid, study, "guide/sense not reverse-complementary")
            continue
        start = tx.find(sense) + 1
        if start < 1:
            exclude(rid, study, "sense not found in accession")
            continue
        overhang = guide[19:]
        core = Core(accession, start, 19, sense, ordered_guide=guide, guide_overhang=overhang)
        try:
            value = float(row["%Inhibition"])
        except (KeyError, ValueError):
            exclude(rid, study, "missing %Inhibition")
            continue
        key = (study, row.get("Gene", "").strip(), accession, core.guide, core.duplex_format)
        records.append({
            "record_id": rid, "study": study, "gene": row.get("Gene", "").strip(),
            "assay": "_".join(x for x in (row.get("Cell", ""), row.get("Concentration", ""),
                                          row.get("Hours", "")) if x) or "unspecified",
            "transcript_id": accession, "target_sha256": sha256(tx), "guide_dna": core.guide,
            "duplex_format": core.duplex_format, "endpoint": "percent_inhibition", "value": value,
            "replicate_group": "", "split": "holdout", "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{row.get('PMID','')}",
            "passenger": "", "guide_overhang": overhang, "mapping_status": "publisher_construct_exact",
            "source_record": row["ID"]})
        candidates.append({"transcript_id": accession, "start": start, "length": 19, "sense": sense,
                           "ordered_guide": guide, "guide_overhang": overhang, "passenger": "",
                           "passenger_overhang": ""})
        by_key[key].append(rid)
        transcripts[accession] = tx
    # Replicate groups only where the exact molecule+condition repeats.
    first = {}
    for key, ids in by_key.items():
        if len(ids) > 1:
            first[key] = ids[0]
    for r, c in zip(records, candidates):
        key = (r["study"], r["gene"], r["transcript_id"], c["ordered_guide"], r["duplex_format"])
        r["replicate_group"] = first.get(key, "")
    return records, candidates, transcripts, exclusions, {"studies": sorted(holdout_studies)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snapshot", type=Path, default=ROOT / "data/raw/sirna_data.xlsx")
    p.add_argument("--fasta", type=Path, default=ROOT / "data/raw/target_gene_sequences.fa")
    p.add_argument("--cohort", type=Path, default=ROOT / "cohort_manifest.json")
    p.add_argument("--out", type=Path, default=ROOT / "data/holdout")
    a = p.parse_args()
    records, candidates, transcripts, exclusions, _ = build(a.snapshot, a.fasta, a.cohort)
    if not records:
        raise ValueError("no holdout records placed")
    out = a.out
    # Prune orphan transcripts so a removed source record cannot leave a stale .fa
    # that run_final then re-pins as if it were still evaluated.
    shutil.rmtree(out / "transcripts", ignore_errors=True)
    (out / "transcripts").mkdir(parents=True, exist_ok=True)
    for accession, tx in sorted(transcripts.items()):
        (out / "transcripts" / f"{safe(accession)}.fa").write_text(f">{accession}\n{tx}\n")
    # Candidate list is unique molecules; repeated measurements stay in measurements.tsv.
    seen, unique_candidates = set(), []
    for c in candidates:
        key = (c["transcript_id"], c["start"], c["length"], c["ordered_guide"], c["guide_overhang"])
        if key in seen:
            continue
        seen.add(key)
        unique_candidates.append(c)
    for name, rows in (("measurements.tsv", records), ("candidates.tsv", unique_candidates)):
        with (out / name).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    (out / "exclusions.json").write_text(json.dumps(exclusions, indent=2) + "\n")
    manifest = {
        "dataset": "weight-independent holdout (publisher snapshot, non-Huesken studies)",
        "split": "holdout", "independent_holdout": True, "endpoint": "percent_inhibition, threshold 70",
        "notes": ["ambiguous reference bases (Lamin A 'N' runs) are replaced in place by 'A', preserving coordinates",
                  "records whose accession is absent from the snapshot FASTA (Takayuki Katoh EGFP, UYB79189) are excluded"],
        "cohort_manifest_sha256": sha256_file(a.cohort),
        "source_hashes": {n: sha256_file(a.snapshot.parent / n)
                          for n in ("sirna_data.xlsx", "target_gene_sequences.fa")},
        "studies": sorted({r["study"] for r in records}),
        "records": len(records), "genes": len({r["gene"] for r in records}),
        "excluded_by_study": dict(collections.Counter(e["study"] for e in exclusions)),
        "transcripts": len(transcripts), "exclusions": len(exclusions),
        "artifacts": {str(rel): sha256_file(out / rel) for rel in
                      ["measurements.tsv", "candidates.tsv", "exclusions.json"]
                      + [f"transcripts/{safe(k)}.fa" for k in sorted(transcripts)]},
    }
    (out / "dataset_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: manifest[k] for k in ("records", "genes", "transcripts", "exclusions")}))


if __name__ == "__main__":
    main()
