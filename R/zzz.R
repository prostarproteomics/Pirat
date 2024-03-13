#' @title xxx
#' @description xxxx
#' # https://community.rstudio.com/t/when-to-use-onload-vs-onattach/21953
#' Usually you want .onLoad, which—as the name suggests—runs when the package is 
#' loaded. If something has to happen before anything is run, that's the way to 
#' go. onAttach only runs when the library is attached, e.g. when somebody calls
#' library(your_package). onLoad will also run when somebody loads but doesn't 
#' attach your package by calling your_package::your_function.
#' @return NA
#' @docType package
#' @aliases Pirat-package
#' @name Pirat
"_PACKAGE"



msg <- paste0("This is Pirat v", utils::packageVersion("Pirat"))
packageStartupMessage(msg)

.onLoad <- function(libname, pkgname) {
  
  envPirat <- basilisk::BasiliskEnvironment("envPirat",
    pkgname = "Pirat",
    packages = c("numpy==1.20.2",
      "matplotlib==3.7",
      "pytorch==1.10.0",
      "python==3.9.5"),
    channels = c("pytorch", "stable", "torch"),
    path = "test_dummy")
}
