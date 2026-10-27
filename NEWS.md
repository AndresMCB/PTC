# PTC 1.1.2

The results of the analyses in the paper (inst/PTC.MainApplication.R, TCGA
BRCA and single-cell data) are exactly the same as with PTC 1.1.1. The test
suite in tests/ checks this.

## Bug fixes

* `PTC.RankByContext()` sorted the scores as text. Now the scores are sorted
  as numbers, from the most negative, and pairs without context score (`NaN`)
  are ranked last. The order for the published analyses does not change,
  because all their scores are negative.
* `PTC()` failed when no miRNA-mRNA interaction was found. Now it returns an
  empty `Summary`.
* When the invariance test fails for a gene (e.g. a singular design matrix),
  `PTC()` now gives a warning with the gene name. Before, the gene was dropped
  silently. In the published single-cell analysis this happens for BZW1,
  EIF4G2 and SUB1.
* In inst/PTC.MainApplication.R, the Wanderlust pseudotime is now passed as a
  vector (`unlist(wandTime)`). In current R versions `order()` of a data frame
  gives an error.

## Dependencies

* CancerSubtypes and seqICP are no longer required. The MAD selection
  (`FSbyMAD`) and the decoupled invariance test (`seqICP.s`) are now
  implemented in PTC (R/getDatabyMAD.R and R/decoupledTest.R). CancerSubtypes
  is deprecated in Bioconductor and seqICP has not been updated since 2017.
  PTC now only depends on miRBaseConverter.
* PTC no longer creates an Rplots.pdf file (the histogram drawn by
  `FSbyMAD`).

## Speed

* `PTC()` is about 2.5 times faster (TCGA BRCA example: 31 s to 12 s on our
  test computer), with the same results:
  * `PTC.RankByContext()` converts miRNA names only for the genes being
    ranked.
  * The bootstrap of the decoupled test does not use the n x n residual
    matrix. The p-values are identical to `seqICP.s`.

## Other changes

* miRBaseConverter no longer prints the table of miRNAs with multiple matches.
* Datasets are compressed better (installed size 81 MB to 55 MB).
* Documentation corrections (`silent` and `TScan` arguments, links, examples,
  typos) and new user manual PTC_1.1.2.pdf.
* `citation("PTC")` gives the reference of the paper in Bioinformatics.
* Added tests (tests/testthat).
