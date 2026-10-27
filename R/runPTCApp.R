#' @title runPTCApp
#'
#' @description Opens the PTC graphical interface in a web browser. The interface lets
#' you run PTC on the example datasets or on your own matched miRNA-mRNA expression
#' data (CSV files), change the parameters, and download the inferred interactions.
#' \cr The interface calls the function \link{PTC}, so results are the same as
#' running \link{PTC} from R with the same data and parameters.
#' @param ... Further arguments passed to \code{shiny::runApp} (e.g. \code{port} or
#' \code{launch.browser}).
#' @details The packages \code{shiny}, \code{bslib} and \code{DT} are required.
#' They can be installed with \code{install.packages(c("shiny", "bslib", "DT"))}.
#' @export
#' @seealso \link{PTC}
#' @examples \dontrun{
#' runPTCApp()
#' }

runPTCApp<-function(...){
  needed<-c("shiny","bslib","DT")
  missing<-needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing)>0){
    stop("The PTC interface needs the packages: ", toString(missing),
         ". Please install them with install.packages(c(\"",
         paste(missing, collapse = "\", \""), "\"))", call. = FALSE)
  }
  shiny::runApp(system.file("app", package = "PTC"), ...)
}
