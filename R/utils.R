#' @title xxxxx
#' @description xxx
#' 
#' @export
#' 
#' @return NA
#' 
source_own_pyScripts <- function(){
    
    tryCatch({
      dir.backup <- getwd()
      setwd(system.file(".", package="Pirat"))
      custom_scripts <- c("LBFGS.py", "llk_maximize.py")
      for (i in custom_scripts)
        reticulate::source_python(system.file("python", i, package = "Pirat"))
      setwd(dir.backup)
      },
      warning = function(w) w,
      error = function(e) e
      )
}





