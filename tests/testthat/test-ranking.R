test_that("interactions are ranked from the most negative context score", {
  rel <- cbind(miR = c("a", "b", "c", "d"), mR = c("A", "B", "C", "D"))
  ctx <- data.frame(miRNA = c("a", "b", "c"), Gene.Symbol = c("A", "B", "C"),
                    context...score = c(-0.05, -0.9, 0.2))
  ranked <- suppressMessages(suppressWarnings(PTC.RankByContext(ctx, rel)))
  expect_equal(unname(ranked[, "miR"]), c("b", "a", "c", "d"))
  expect_equal(unname(ranked[, "Score"]), c("-0.9", "-0.05", "0.2", "NaN"))
})

test_that("no interactions gives an empty Summary instead of an error", {
  empty <- InterList.toMatrix(list())
  expect_equal(nrow(empty), 0)
  ranked <- PTC.RankByContext(TS7.0_Conserved_Site_Context_Scores, empty)
  expect_equal(nrow(ranked), 0)
  expect_equal(colnames(ranked), c("rank", "miR", "mR", "Score"))
})
