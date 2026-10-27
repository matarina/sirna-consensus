import random
import sys
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

TRANSCRIPT = (
    "CACCCTACAATTTAATCTCGTCCTGTCGGGAAACGGCGAGCTGGATTGGTGCACTGGCAACATCCCGCGCTAGGTTTGATGGTAAACTGTACGAGCCACC"
    "CTGCTGCAGCGTCCACGGACTTACAGTAAGCGCCCTCGAGGCTCCTCCCTCCATTGAAAGCAGAAGTCCCGGATTTCTTAACAATCTACCTCTGGCACGC"
)
TRANSCRIPT_ID = "TEST_0001"


@pytest.fixture
def transcript() -> str:
    return TRANSCRIPT


@pytest.fixture
def transcript_id() -> str:
    return TRANSCRIPT_ID


@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
    path = tmp_path / "test.fa"
    path.write_text(f">{TRANSCRIPT_ID} synthetic\n{TRANSCRIPT}\n")
    return path


@pytest.fixture
def random_dna():
    def make(length: int, seed: int = 7) -> str:
        rng = random.Random(seed)
        return "".join(rng.choice("ACGT") for _ in range(length))
    return make
