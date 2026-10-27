# PTC graphical interface. Start it with PTC::runPTCApp().
# The interface prepares the data and calls PTC(), so results are the same as
# running PTC() from R with the same data and parameters.

library(shiny)
library(bslib)
library(PTC)

options(shiny.maxRequestSize = 500 * 1024^2)  # uploads up to 500 MB

# Reads a CSV with sample names in the first column and gene names in the header.
# If genes.in.rows is TRUE the table is transposed (samples in columns).
readExpression <- function(path, genes.in.rows) {
  d <- read.csv(path, row.names = 1, check.names = FALSE)
  d <- as.matrix(d)
  if (genes.in.rows) d <- t(d)
  if (!is.numeric(d)) {
    stop("The expression values must be numbers. Please check the file format.")
  }
  d
}

# Reads a CSV with two columns: sample name and pseudotime value.
readPseudotime <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  if (ncol(d) < 2) {
    stop("The pseudotime file must have two columns: sample name and pseudotime.")
  }
  setNames(as.numeric(d[[2]]), as.character(d[[1]]))
}

ui <- page_sidebar(
  title = "PTC: Pseudo-temporal Causality",
  theme = bs_theme(version = 5, primary = "#2c5d8f"),
  sidebar = sidebar(
    width = 360,
    h5("1. Data"),
    radioButtons("source", NULL,
                 c("Example: breast cancer bulk data (TCGA BRCA)" = "tcga",
                   "Example: single-cell data (Wang et al., 2019)" = "sc",
                   "My own data" = "upload")),
    conditionalPanel(
      "input.source == 'upload'",
      fileInput("mirFile", "miRNA expression (CSV)", accept = ".csv"),
      fileInput("mrFile", "mRNA expression (CSV)", accept = ".csv"),
      checkboxInput("genesInRows", "Genes are in rows (samples in columns)", FALSE),
      helpText("First column: sample names. Header: gene names.",
               "Samples are matched by name, or by row order if the names differ.")
    ),
    h5("2. Pseudotime"),
    radioButtons("timeSource", NULL,
                 c("Order samples by a marker gene" = "marker",
                   "Upload a pseudotime" = "file")),
    conditionalPanel(
      "input.timeSource == 'marker'",
      textInput("marker", "Marker gene", "VIM")
    ),
    conditionalPanel(
      "input.timeSource == 'file'",
      fileInput("timeFile", "Pseudotime (CSV: sample, value)", accept = ".csv")
    ),
    h5("3. Parameters"),
    layout_columns(
      numericInput("nmiR", "miRNAs", 30, min = 1, step = 1),
      numericInput("nmR", "mRNAs", 1500, min = 1, step = 1)
    ),
    layout_columns(
      numericInput("alpha", "alpha", 0.02, min = 0.0001, max = 0.5, step = 0.01),
      numericInput("ngrid", "Segments", 2, min = 2, step = 1)
    ),
    checkboxInput("complements", "Compare each segment against its complement", TRUE),
    checkboxInput("exploreAll", "Explore all sets of miRNAs", TRUE),
    accordion(
      open = FALSE,
      accordion_panel(
        "Advanced: own TargetScan predictions",
        fileInput("tsFile", "miRNA-mRNA predictions (CSV)", accept = ".csv"),
        helpText("Two columns: miRNA (miRBase v21 names, e.g. hsa-miR-21-5p)",
                 "and mRNA (gene symbol). If empty, TargetScan 7.0 is used.")
      )
    ),
    actionButton("run", "Run PTC", class = "btn-primary btn-lg")
  ),

  uiOutput("main")
)

server <- function(input, output, session) {
  data("GroundT", package = "PTC", envir = environment())

  # Expression data of the selected source
  expression <- reactive({
    if (input$source == "tcga") {
      data("TCGA_BRCAdata", package = "PTC", envir = environment())
      list(miRs = TCGA_BRCAdata$miRs, mRNAs = TCGA_BRCAdata$mRNAs)
    } else if (input$source == "sc") {
      PTC:::prepareSCdata()
    } else {
      validate(need(input$mirFile, "Please upload the miRNA expression file."),
               need(input$mrFile, "Please upload the mRNA expression file."))
      miRs <- readExpression(input$mirFile$datapath, input$genesInRows)
      mRNAs <- readExpression(input$mrFile$datapath, input$genesInRows)
      # Match samples by name; if the names differ, by row order (as PTC does)
      samples <- intersect(rownames(miRs), rownames(mRNAs))
      if (length(samples) >= 2 * input$ngrid) {
        dropped <- length(union(rownames(miRs), rownames(mRNAs))) - length(samples)
        note <- if (dropped > 0) {
          paste(length(samples), "samples are in both files and were used;",
                dropped, "samples in only one file were left out.")
        }
        list(miRs = miRs[samples, , drop = FALSE], mRNAs = mRNAs[samples, , drop = FALSE],
             byName = TRUE, note = note)
      } else {
        validate(need(nrow(miRs) == nrow(mRNAs),
                      paste("The miRNA and mRNA files have different samples.",
                            "Use the same sample names in both files, or list the",
                            "samples in the same order.")))
        list(miRs = miRs, mRNAs = mRNAs, byName = FALSE,
             note = paste("The sample names of the two files are different, so samples",
                          "were matched by row order: row 1 of the miRNA file with row 1",
                          "of the mRNA file, and so on."))
      }
    }
  })

  result <- reactiveVal(NULL)
  failed <- reactiveVal(character(0))

  observeEvent(input$run, {
    result(NULL)
    failed(character(0))
    d <- tryCatch(expression(), error = function(e) {
      showNotification(conditionMessage(e), type = "error", duration = NULL)
      NULL
    })
    if (is.null(d)) return()

    # Pseudotime
    if (input$timeSource == "marker") {
      marker <- toupper(trimws(input$marker))
      genes <- toupper(colnames(d$mRNAs))
      if (!marker %in% genes) {
        showNotification(paste0("The marker gene ", marker, " is not in the mRNA data."),
                         type = "error", duration = NULL)
        return()
      }
      VIM <- d$mRNAs[, match(marker, genes)]
    } else {
      if (is.null(input$timeFile)) {
        showNotification("Please upload the pseudotime file.", type = "error")
        return()
      }
      pt <- readPseudotime(input$timeFile$datapath)
      if (!all(rownames(d$mRNAs) %in% names(pt))) {
        showNotification("The pseudotime file must contain a value for every sample.",
                         type = "error", duration = NULL)
        return()
      }
      VIM <- pt[rownames(d$mRNAs)]
    }

    TScan <- NULL
    if (!is.null(input$tsFile)) {
      TScan <- as.matrix(read.csv(input$tsFile$datapath, check.names = FALSE))[, 1:2]
    }

    nmiR <- min(input$nmiR, ncol(d$miRs))
    nmR <- min(input$nmR, ncol(d$mRNAs))
    warnings <- character(0)

    res <- withProgress(message = "Running PTC", value = 0, {
      tryCatch(
        withCallingHandlers(
          PTC(miRNAs = d$miRs, mRNAs = d$mRNAs, VIM = VIM,
              nmiR = nmiR, nmR = nmR, ngrid = input$ngrid, alpha = input$alpha,
              complements = input$complements, explore.all = input$exploreAll,
              TScan = TScan),
          message = function(m) {
            gene <- sub(".*Current mRNA  = \\{ (.*) \\}.*", "\\1", conditionMessage(m))
            if (gene != conditionMessage(m)) incProgress(1 / nmR, detail = gene)
            invokeRestart("muffleMessage")
          },
          warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
          }),
        error = function(e) {
          showNotification(paste("PTC stopped:", conditionMessage(e)),
                           type = "error", duration = NULL)
          NULL
        })
    })
    if (is.null(res)) return()

    res$note <- d$note
    res$settings <- list(source = input$source, byName = !isFALSE(d$byName), timeSource = input$timeSource,
                         marker = input$marker, nmiR = nmiR, nmR = nmR,
                         ngrid = input$ngrid, alpha = input$alpha,
                         complements = input$complements,
                         exploreAll = input$exploreAll, ownTScan = !is.null(TScan))
    failed(grep("invariance test failed", warnings, value = TRUE))
    result(res)
  })

  # Results table: PTC Summary plus whether the interaction is experimentally confirmed
  resultTable <- reactive({
    res <- result()
    req(res)
    s <- res$Summary
    data.frame(Rank = as.integer(s[, "rank"]),
               miRNA = s[, "miR"],
               mRNA = s[, "mR"],
               Score = as.numeric(s[, "Score"]),
               Confirmed = ifelse(mapply(function(mi, mr) mi %in% GroundT[[mr]],
                                         s[, "miR"], s[, "mR"]), "yes", "no"),
               row.names = NULL)
  })

  output$main <- renderUI({
    if (is.null(result())) {
      return(card(
        card_header("Pseudo-temporal Causality (PTC)"),
        p("PTC finds causal miRNA-mRNA regulatory relationships during a biological",
          "process. Samples are ordered by pseudotime (by default, the expression of",
          "the EMT marker VIM), and a miRNA is reported as a cause of a target mRNA when",
          "the linear model between them stays invariant along the pseudotime."),
        tags$ol(
          tags$li("Choose an example dataset or upload your matched miRNA and mRNA expression."),
          tags$li("Choose how to order the samples (pseudotime)."),
          tags$li("Keep the default parameters (used in the paper) or change them."),
          tags$li("Click ", tags$b("Run PTC"), ". The breast cancer example takes about",
                  "half a minute and the single-cell example a few seconds.")
        ),
        p(tags$b("Citation: "), "Cifuentes-Bernal AM, Pham VVH, Li X, Liu L, Li J, Le TD.",
          "A pseudo-temporal causality approach to identifying miRNA-mRNA interactions",
          "during biological processes. Bioinformatics, 2021.",
          tags$a("doi:10.1093/bioinformatics/btaa899",
                 href = "https://doi.org/10.1093/bioinformatics/btaa899", target = "_blank"))
      ))
    }
    tbl <- resultTable()
    tagList(
      layout_columns(
        value_box("Interactions", nrow(tbl)),
        value_box("Target mRNAs", length(unique(tbl$mRNA))),
        value_box("Experimentally confirmed", sum(tbl$Confirmed == "yes"),
                  p("in miRTarBase, TarBase or miRWalk"))
      ),
      if (!is.null(result()$note)) {
        card(class = "border-info", card_body(result()$note))
      },
      if (length(failed()) > 0) {
        card(class = "border-warning",
             card_header("The invariance test could not be computed for some genes"),
             p("These genes are reported without parents:"),
             tags$ul(lapply(sub("The invariance test failed for (.*?); .*", "\\1", failed()),
                            tags$li)))
      },
      card(
        card_header(
          class = "d-flex justify-content-between align-items-center",
          "Inferred miRNA-mRNA interactions",
          downloadButton("download", "Download CSV", class = "btn-sm")
        ),
        DT::DTOutput("table")
      ),
      card(
        card_header("Reproduce this analysis in R"),
        tags$pre(rCode())
      )
    )
  })

  output$table <- DT::renderDT({
    DT::datatable(resultTable(), rownames = FALSE,
                  options = list(pageLength = 25, order = list(list(0, "asc"))))
  })

  output$download <- downloadHandler(
    filename = function() paste0("PTC_interactions_", Sys.Date(), ".csv"),
    content = function(file) write.csv(resultTable(), file, row.names = FALSE)
  )

  # R code that gives the same result with the package
  rCode <- reactive({
    s <- result()$settings
    data <- switch(s$source,
      tcga = c("data(TCGA_BRCAdata)",
               "miRNAs <- TCGA_BRCAdata$miRs",
               "mRNAs <- TCGA_BRCAdata$mRNAs"),
      sc = c("SC <- PTC:::prepareSCdata()  # as in inst/PTC.MainApplication.R",
             "miRNAs <- SC$miRs",
             "mRNAs <- SC$mRNAs"),
      upload = c(
        sprintf('miRNAs <- as.matrix(read.csv("%s", row.names = 1, check.names = FALSE))',
                input$mirFile$name),
        sprintf('mRNAs <- as.matrix(read.csv("%s", row.names = 1, check.names = FALSE))',
                input$mrFile$name),
        if (input$genesInRows) c("miRNAs <- t(miRNAs)", "mRNAs <- t(mRNAs)"),
        if (s$byName) c("samples <- intersect(rownames(miRNAs), rownames(mRNAs))",
                        "miRNAs <- miRNAs[samples, ]",
                        "mRNAs <- mRNAs[samples, ]")))
    time <- if (s$timeSource == "marker") {
      sprintf('pseudotime <- mRNAs[, "%s"]', toupper(trimws(s$marker)))
    } else {
      c(sprintf('pt <- read.csv("%s")', input$timeFile$name),
        "pseudotime <- setNames(pt[[2]], pt[[1]])[rownames(mRNAs)]")
    }
    paste(c("library(PTC)", data, time,
            if (s$ownTScan) sprintf('TScan <- as.matrix(read.csv("%s"))[, 1:2]', input$tsFile$name),
            sprintf(paste0("result <- PTC(miRNAs, mRNAs, VIM = pseudotime, nmiR = %d, nmR = %d,\n",
                           "              ngrid = %d, alpha = %s, complements = %s,\n",
                           "              explore.all = %s%s)"),
                    as.integer(s$nmiR), as.integer(s$nmR), as.integer(s$ngrid), s$alpha,
                    s$complements, s$exploreAll, if (s$ownTScan) ", TScan = TScan" else ""),
            "result$Summary"),
          collapse = "\n")
  })
}

shinyApp(ui, server)
