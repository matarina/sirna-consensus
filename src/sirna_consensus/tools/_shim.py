"""Run a legacy single-file script under a compatibility shim.

``siRNADesigner.py`` was written for a Biopython where ``Seq`` behaved like
``str``. On modern Biopython two things break:

1. ``re`` refuses ``Seq`` objects (it needs a real ``str`` subclass);
2. ``SeqIO.parse`` hands back ``Seq`` records that then get sliced and joined.

The script also never closes its own output file, so under ``runpy`` the buffer
is only flushed by luck. We do all three here and leave upstream untouched.

Usage:  python _shim.py /path/to/script.py <script args...>
"""
from __future__ import annotations

import gc
import io
import re
import runpy
import sys
from typing import Any


def _patch_seqio() -> None:
    from Bio import SeqIO

    class Record:
        __slots__ = ("seq", "id", "name", "description")

        def __init__(self, record: Any) -> None:
            self.seq = str(record.seq)
            self.id = str(record.id)
            self.name = str(getattr(record, "name", record.id))
            self.description = str(getattr(record, "description", ""))

    original = SeqIO.parse

    def parse(*args: Any, **kwargs: Any):
        for record in original(*args, **kwargs):
            yield Record(record)

    SeqIO.parse = parse


def _patch_re() -> None:
    for name in ("match", "search", "fullmatch", "findall", "finditer", "split"):
        original = getattr(re, name)

        def make(original):
            def wrapper(pattern, string, *args, **kwargs):
                if string is not None and not isinstance(string, (str, bytes)):
                    string = str(string)
                return original(pattern, string, *args, **kwargs)
            return wrapper

        setattr(re, name, make(original))


def main() -> None:
    _patch_seqio()
    _patch_re()
    script = sys.argv[1]
    sys.argv = sys.argv[1:]
    namespace = runpy.run_path(script, run_name="__main__")
    for value in list(namespace.values()):
        if isinstance(value, io.IOBase):
            try:
                value.close()
            except Exception:  # noqa: BLE001 - best effort flush
                pass
    gc.collect()


if __name__ == "__main__":
    main()
