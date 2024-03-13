#' @title xxx
#' @rdname install_pirat
#' 
#' @examples
#' install_pirat()
#' @export
#' 
#' @importFrom stats cov
#' 
#' @return NA
#' 
install_pirat <- function() { 
  proc <- basilisk::basiliskStart(envPirat)
  on.exit(basilisk::basiliskStop(proc))
  
  some_useful_thing <- basilisk::basiliskRun(proc, 
    fun = function() {
      py <- reticulate::import("torch", delay_load = FALSE)
      
      message('Installation completed !')
      #output <- NULL
      # The return value MUST be a pure R object, i.e., no reticulate
      # Python objects, no pointers to shared memory. 
      #output 
    })
  
  some_useful_thing
}