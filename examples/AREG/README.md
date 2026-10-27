# AREG example run

Real run against NM_001657.4 (1234 nt) with all three tools installed.
Reproduce with:

```bash
docker pull ghcr.io/austin-s-h/sirnaforge:sha-4c05e4b   # or ghcr.nju.edu.cn mirror
sirna-consensus run --fasta NM_001657.4.fa --out examples/AREG \
  --tools sirnadesigner,siren,sirnaforge --exclude-singletons
```

Result: 1246 hits placed from 3 tools, 898 distinct sites.
Of those, 29 were proposed by all three tools (3/3) and 290 by two (2/3).
174 sites passed the sequence filters; the 10 selected span >=30 nt.

Note the vote distribution: unanimous sites are a small fraction. That is the
expected behaviour, and the reason the full site list is kept in
candidates.tsv rather than only the intersection.
