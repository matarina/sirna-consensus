# Sources and exposure audit

- Huesken et al., 2005, DOI 10.1038/nbt1118: original publisher supplementary
  tables and source hashes already curated in the sibling project. Reused as
  exposed development data. Original constructs are the coordinate reference.
- siRNAEfficacyDB publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC11665841/
- Publisher-linked Figshare dataset: https://doi.org/10.6084/m9.figshare.25908841.v1
  Files 46578241 and 46578250. Snapshot contains 3,544 measurements from seven
  source studies. Reproduce acquisition with `scripts/fetch_snapshot.py`.
- OligoFormer source: https://github.com/lulab/OligoFormer . Local model and
  thermodynamic sources/checkpoints are pinned in the scorer manifest. The seven
  database source studies appear in the published development/evaluation
  collections; released checkpoint-specific splits are not verified.
- RNA-FM: https://github.com/ml4bio/RNA-FM and
  https://huggingface.co/cuhkaih/rnafm . Runtime uses RNA-FM 0.2.2 and a pinned local
  `RNA-FM_pretrained.pth` checkpoint, with no runtime network downloads.
- Specificity reference: existing Ensembl release 116 GRCh38 cDNA and Bowtie 1.3.1
  indexes, validated against the sibling reference manifest.

The Figshare metadata says CC BY 4.0, whereas the database paper describes
CC BY-NC 2.0 UK. Preserve both statements and keep raw snapshots local pending
reconciliation. Pipeline code remains MIT; external models retain their licenses.

The snapshot's Huesken cell/assay metadata conflicts with original-source
curation. Its `%Inhibition` column rescales normalized reporter activity. The
benchmark retains original H1299/fluorescent-reporter metadata and raw normalized
activity, and does not silently convert these data to endogenous percent knockdown.

No independent cohort is currently established. Neither changing a random seed
nor repartitioning the already exposed Huesken genes establishes one. Independent
source/checkpoint exposure review is a prerequisite for final validation.
