test_that("the interface gives the same results as PTC()", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("bslib")
  skip_if_not_installed("DT")
  expected <- readRDS(test_path("published_sc.rds"))$Summary
  app <- shiny::shinyAppDir(system.file("app", package = "PTC"))

  shiny::testServer(app, {
    session$setInputs(source = "sc", timeSource = "marker", marker = "vim",
                      nmiR = 30, nmR = 1500, alpha = 0.02, ngrid = 2,
                      complements = TRUE, exploreAll = TRUE, genesInRows = FALSE)
    session$setInputs(run = 1)
    expect_identical(result()$Summary, expected)
    expect_equal(nrow(resultTable()), nrow(result()$Summary))
    expect_length(failed(), 3)  # BZW1, EIF4G2 and SUB1, as in the paper
  })
})
