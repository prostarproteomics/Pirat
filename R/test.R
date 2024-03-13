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


my_example_function <- function(ARG_VALUE_1, ARG_VALUE_2) { 
    proc <- basiliskStart(my_env)
    source_own_pyScripts()
    on.exit(basiliskStop(proc))

    some_useful_thing <- basiliskRun(proc, fun=function(arg1, arg2) {
        mod <- reticulate::import("torch")
        output <- 3

        # The return value MUST be a pure R object, i.e., no reticulate
        # Python objects, no pointers to shared memory. 
        output 
    }, arg1=ARG_VALUE_1, arg2=ARG_VALUE_2)

    some_useful_thing
}
