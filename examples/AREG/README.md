# AREG example run

Real run against NM_001657.4 (1234 nt) with two of the three tools installed
(sirnaforge needs its OCI image; see README). Reproduce with:

```bash
sirna-consensus run --fasta NM_001657.4.fa --out examples/AREG \
  --tools sirnadesigner,siren --exclude-singletons
```

Result: 687 hits placed, 652 sites, 35 unanimous (2/2), 16 passing the
sequence filters, 10 selected with >=30 nt spacing.
