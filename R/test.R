my_example_function <- function(ARG_VALUE_1, ARG_VALUE_2) { 
    proc <- basiliskStart(my_env)
    on.exit(basiliskStop(proc))

    some_useful_thing <- basiliskRun(proc, fun=function(arg1, arg2) {
        mod <- reticulate::import("scikit-learn")
        output <- mod$some_calculation(arg1, arg2)

        # The return value MUST be a pure R object, i.e., no reticulate
        # Python objects, no pointers to shared memory. 
        output 
    }, arg1=ARG_VALUE_1, arg2=ARG_VALUE_2)

    some_useful_thing
}
