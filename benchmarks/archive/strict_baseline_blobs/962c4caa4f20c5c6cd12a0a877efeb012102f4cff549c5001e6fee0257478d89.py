"""verify_freezes must report a bad lock, not crash, and count missing entries."""
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import verify_freezes  # noqa: E402


def test_malformed_lock_raises_valueerror(tmp_path):
    lock = tmp_path / "l.json"
    lock.write_text("{ not json")
    with pytest.raises(ValueError):
        verify_freezes.check(lock)


def test_non_dict_artifacts_raises_valueerror(tmp_path):
    lock = tmp_path / "l.json"
    lock.write_text("[]")                       # json.loads -> list, ["artifacts"] is TypeError
    with pytest.raises(ValueError):
        verify_freezes.check(lock)
    lock.write_text(json.dumps({"artifacts": "nope"}))
    with pytest.raises(ValueError):
        verify_freezes.check(lock)


def test_missing_and_directory_are_counted_not_crashed(tmp_path):
    (tmp_path / "sub").mkdir()
    (tmp_path / "a.txt").write_text("x")
    lock = tmp_path / "l.json"
    lock.write_text(json.dumps({"artifacts": {"a.txt": "deadbeef", "sub": "x", "gone": "x"}}))
    present, mismatched, missing = verify_freezes.check(lock)
    assert present == 0 and mismatched == ["a.txt"] and missing == 2   # dir + absent


def test_main_reports_bad_lock_and_continues(tmp_path, capsys):
    bad = tmp_path / "bad.json"
    bad.write_text("{ not json")
    good = tmp_path / "good.json"
    (tmp_path / "ok.txt").write_text("hi")
    import hashlib
    good.write_text(json.dumps({"artifacts": {"ok.txt": hashlib.sha256(b"hi").hexdigest()}}))
    assert verify_freezes.main([str(bad), str(good)]) == 1   # failed on bad, still checked good
    out = capsys.readouterr().out
    assert "NOT PARSEABLE" in out and "good.json: ok" in out


def test_historical_source_snapshot_is_required_and_checked(tmp_path):
    import hashlib

    source = tmp_path/'source.py'
    source.write_text('historical implementation\n')
    lock = tmp_path/'lock.json'
    key = '../src/module.py'
    lock.write_text(json.dumps({'artifacts':{key:hashlib.sha256(source.read_bytes()).hexdigest()}}))
    mapping = {key:str(source)}
    assert verify_freezes.check(lock,mapping)==(1,[],0)
    source.write_text('modified\n')
    assert verify_freezes.check(lock,mapping)[1]==[key]
    source.unlink()
    assert 'required historical source snapshot missing' in verify_freezes.check(lock,mapping)[1][0]
