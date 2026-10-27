# The expected results (published_*.rds) were produced with PTC 1.1.1, the
# version used for the paper, with the settings of inst/PTC.MainApplication.R.
# Any change to the package must give exactly the same results.

published_parts <- function(res) res[c("Index", "Names", "genes", "Summary")]

test_that("single-cell results are the published ones", {
  SC <- PTC:::prepareSCdata()
  # As in the published analysis, the invariance test cannot be computed for
  # BZW1, EIF4G2 and SUB1 (singular design matrix), so PTC warns and they
  # have no parents.
  res <- suppressWarnings(suppressMessages(
    PTC(miRNAs = SC$miRs, mRNAs = SC$mRNAs, VIM = SC$mRNAs[, "VIM"])))
  expect_identical(published_parts(res),
                   readRDS(test_path("published_sc.rds")))
})

test_that("TCGA BRCA results are the published ones", {
  skip_on_cran()  # takes some seconds
  res <- suppressMessages(PTC(miRNAs = TCGA_BRCAdata$miRs,
                              mRNAs = TCGA_BRCAdata$mRNAs,
                              VIM = TCGA_BRCAdata$mRNAs[, "VIM"]))
  expect_identical(published_parts(res),
                   readRDS(test_path("published_tcga.rds")))
})
