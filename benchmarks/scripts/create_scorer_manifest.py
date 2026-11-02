"""Pin the local frozen scorer and dependency source files.

The manifest never asserts training independence unconditionally: the training
audit text is derived from ``cohort_manifest.json`` when that audit has actually
verified the checkpoint's training exposure.
"""
import argparse
import hashlib
import importlib.metadata
import json
import sys
from pathlib import Path

try:
    import fm
except ImportError:  # the repo interpreter may lack RNA-FM; hashing still works
    fm = None

DOWNSTREAM_SOURCE_CMD = "run benchmarks/scripts/audit_exposure.py"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def training_audit_text(cohort_path):
    if not cohort_path.is_file():
        return "unverified checkpoint split; " + DOWNSTREAM_SOURCE_CMD
    cohort = json.loads(cohort_path.read_text())
    if (cohort.get("checkpoint_training_audit") == "verified"
            and cohort.get("independent_holdout_eligible")
            and cohort.get("holdout_studies")):
        return (f"checkpoint training source exposure audited (upstream Hu.csv); "
                f"{len(cohort['holdout_studies'])} weight-independent holdout studies / "
                f"{cohort['holdout_records']} records / {cohort['holdout_genes']} genes exit every "
                f"training and development molecule; upstream evaluation exposure disclosed. "
                f"See benchmarks/cohort_manifest.json")
    return "cohort unverified: " + (cohort.get("reason") or "no reason recorded")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--upstream", type=Path, required=True)
    p.add_argument("--rna-fm-checkpoint", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--device", default="cpu")
    a = p.parse_args()

    root = Path(__file__).resolve().parents[2]
    worker = root / "benchmarks/scripts/oligoformer_worker.py"
    model = a.upstream / "scripts/model.py"
    features = a.upstream / "scripts/infer.py"
    checkpoint = a.upstream / "model/best_model.pth"
    if fm is None:
        print("warning: RNA-FM not importable here; fm artifacts dropped and versions may be null "
              "(use the scorer interpreter for the canonical manifest)", file=sys.stderr)
    artifacts = [worker, model, features, checkpoint, a.rna_fm_checkpoint]
    artifacts += list(Path(fm.__file__).parent.rglob("*.py")) if fm else []

    def version(name):
        try:
            return importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            return None


    d = {
        "name": "OligoFormer", "input_convention": "paired_antisense_19",
        "supported_formats": ["binding19-guide19-overhang0-passenger_unknown",
                              "binding19-guide21-overhang2-passenger_unknown",
                              "binding19-guide21-overhang2-passenger21-overhang2"],
        "checkpoint_sha256": digest(checkpoint),
        "artifacts": {str(p.resolve()): digest(p) for p in artifacts},
        "command": [sys.executable, str(worker), "--input", "{input}", "--output", "{output}",
                    "--model-source", str(model.resolve()), "--features-source", str(features.resolve()),
                    "--checkpoint", str(checkpoint.resolve()), "--rna-fm-checkpoint",
                    str(a.rna_fm_checkpoint.resolve()), "--device", a.device, "--threads", "4"],
        "packages": {name: version(name) for name in ("torch", "numpy", "pandas", "rna-fm")},
        "inference_policy": "eval mode, fixed seed 2026, batch size one, no fitting, local checkpoint only",
        "training_audit": training_audit_text(root / "benchmarks/cohort_manifest.json"),
        "timeout_s": 3600, "independent_validation_passed": False,
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(d, indent=2) + "\n")
    print(a.out)


if __name__ == "__main__":
    main()
