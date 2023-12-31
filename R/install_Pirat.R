#' @title Install Pirat package
#' 
#' @description This script installs Python and PyTorch in the requested
#' versions. It is largely inspired by wthe scripts in rTorch package 
#' (https://github.com/f0nzie/rTorch)
#'
#' @param method Installation method. By default, "auto" automatically finds a
#'   method that will work in the local environment. Change the default to force
#'   a specific installation method.  Note that since this command runs without 
#'   privilege the "system" method is available only on _Windows_.
#'
#' @param conda www
#' 
#' @param version PyTorch version to install. The "default" version is 
#' __1.10.0__. You can specify a specific __PyTorch__ version with 
#' `version="1.2"`, or `version="1.6"`.
#'
#' @param envname Name of Python or conda environment to install within.
#'   The default environment name is `r-pirat`.
#'
#' @param extra_packages Additional Python packages to install along with
#'   PyTorch. Default are `c("numpy=1.20.2", "matplotlib")`.
#'
#' @param restart_session Restart R session after installing (note this will
#'   only occur within RStudio).
#'
#' @param conda_python_version the _Python_ version installed in the created _conda_
#'   environment. Python __3.9.5__ is installed by default. But you could 
#'   specify for instance: `conda_python_version="3.7"`.
#'
#' @param pip logical
#'
#' @param channel conda channel. The default channel is `stable`.
#'   The alternative channel is `nightly`.
#'
#' @param cuda_version string for the cuda toolkit version to install. For example,
#'   to install a specific CUDA version use `cuda_version="10.2"`.
#'
#' @param dry_run logical, set to TRUE for unit tests, otherwise will execute
#'   the command.
#'
#' @param ... other arguments passed to [reticulate::conda_install()] or
#'   [reticulate::virtualenv_install()].
#'
#' @importFrom jsonlite fromJSON
#' @examples
#' \dontrun{
#'
#' # install PyTorch 1.10.0 on Python 3.3.9.5 including pandas
#' install_pytorch(version = "1.10.0", conda_python_version = "3.9.5",
#' extra_packages = "pandas")
#'
#' # Install PyTorch 1.10.0, Python 3.9.5, pandas, matplotlib install from the console
#' install_pytorch(version = "1.10.0", conda_python_version = "3.9.5",
#' extra_packages = c("pandas", "matplotlib"))
#'
#' # Install PyTorch 1.10.0 on Python 3.9.5 including pandas, matplotlib
#' install_pytorch(version = "1.10.0", conda_python_version = "3.9.5",
#' extra_packages = c("pandas", "matplotlib"), dry_run = FALSE)
#' }
#'
#' @export
#' 
install_pirat <- function(method = c("conda", "virtualenv", "auto"),
                          conda = "auto",
                          version = requested_versions$torch,
                          envname = pirat_envname,
                          extra_packages = NULL,
                          restart_session = TRUE,
                          conda_python_version = requested_versions$python,
                          pip = FALSE,
                          channel = "stable",
                          cuda_version = NULL,
                          dry_run = FALSE,
                          ...) {
  
  
  # verify 64-bit
  if (.Machine$sizeof.pointer != 8) {
    stop("Unable to install PyTorch on this platform.",
         "Binary installation is only available for 64-bit platforms.")
  }
  
  method <- match.arg(method)
  
  # unroll version
  ver <- parse_torch_version(version, cuda_version, channel)
  
  version <- ver$version
  gpu <- ver$gpu
  package <- ver$package
  cpu_gpu_packages <- ver$cpu_gpu_packages
  channel <- ver$channel
  
  # Packages in this list should always be installed.
  
  default_packages <- c("numpy=1.20.2", 'matplotlib')
  
  # # Resolve torch probability version.
  # if (!is.na(version) && substr(version, 1, 4) %in% c("1.1.0", "1.1", "1.1.0")) {
  #   default_packages <- c(default_packages, "pandas")
  #   # install pytorch-nightly
  # } else if (is.na(version) ||(substr(version, 1, 4) %in% c("2.0.") || version == "nightly")) {
  #   default_packages <- c(default_packages, "numpy")
  # }
  
  extra_packages <- unique(c(cpu_gpu_packages, default_packages, extra_packages))
  
  if (dry_run) {
    os <- ifelse(is_osx(), "osx",
                 ifelse(is_linux(), "linux",
                        ifelse(is_windows(), "windows", "None")))
    out <- list(package = package, 
                extra_packages = extra_packages,
                envname = envname, 
                conda = conda,
                conda_python_version = 
                conda_python_version,
                channel = channel, 
                pip = pip, 
                os = os)
    return(out)
  }
  
  # Main OS verification.
  if (is_osx() || is_linux()) {
    
    if (method == "conda") {
      install_conda(
        package = package,
        extra_packages = extra_packages,
        envname = envname,
        conda = conda,
        conda_python_version = conda_python_version,
        channel = channel,
        pip = pip,
        ...
      )
    } else if (method == "virtualenv" || method == "auto") {
      install_virtualenv(
        package = package,
        extra_packages = extra_packages,
        envname = envname,
        ...
      )
    }
    
  } else if (is_windows()) {
    
    if (method == "virtualenv") {
      stop("Installing PyTorch into a virtualenv is not supported on Windows",
           call. = FALSE)
    } else if (method == "conda" || method == "auto") {
      
      install_conda(
        package = package,
        extra_packages = extra_packages,
        envname = envname,
        conda = conda,
        conda_python_version = conda_python_version,
        channel = channel,
        pip = pip,
        ...
      )
      
    }
    
  } else {
    stop("Unable to install PyTorch on this platform. ",
         "Binary installation is available for Windows, OS X, and Linux")
  }
  
  message("\nInstallation complete.\n\n")
  
  
  is.rstudio <- function(){
    .Platform$GUI == "RStudio"
  }
  
  if (restart_session){
    if (is.rstudio() &&
        requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::hasFun("restartSession"))
      rstudioapi::restartSession(command='library(Pirat)')
    else
      cat("Please restart the R session and reload the 'Pirat' package.")
  }
  
  invisible(NULL)
}



install_conda <- function(package, 
                          extra_packages, 
                          envname, 
                          conda,
                          conda_python_version, 
                          channel, 
                          pip, 
                          ...) {
  

  remove_Pirat(envname, conda)
  
  
  message("Creating ", envname, " conda environment... \n")
  reticulate::conda_create(
    envname = envname, conda = conda,
    packages = paste0("python=", conda_python_version)
  )
  
  message("Installing python modules...\n")
  # rTorch::conda_install(envname="r-torch-37", packages="pytorch-cpu",
  #         channel = "pytorch", conda="auto", python_version = "3.7")
  conda_install(
    envname = envname,
    packages = c(package, extra_packages),
    conda = conda,
    pip = pip,       # always use pip since it's the recommend way.
    channel = channel,
    ...
  )
  
}

install_virtualenv <- function(package, extra_packages, envname, ...) {
  
  remove_Pirat(envname)
  
  
  message("Creating ", envname, " virtualenv environment... \n")
  reticulate::virtualenv_create(envname = envname)
  
  message("Installing python modules...\n")
  reticulate::virtualenv_install(
    envname = envname,
    packages = c(package, extra_packages),
    ...
  )
  
}
  



parse_torch_version <- function(version, cuda_version = NULL, channel = "stable") {
  default_version <- "1.10.0"
  # channel <- "pytorch"    # this is the channel
  
  ver <- list(
    version = default_version,
    gpu = FALSE,
    package = NULL,
    cuda_version = cuda_version,
    cpu_gpu_packages = NULL,
    channel = channel
  )
  
  if (version == "default") {
    ver$package <- paste0("pytorch==", ver$version)
  } else {
    ver$version <- version
    ver$package <- paste0("pytorch==", ver$version)
  }
  
  
  if (is.null(ver$cuda_version)) {
    ver$cpu_gpu_packages <- "cpuonly"
  } else {
    ver$cuda_version <- cuda_version
    ver$cpu_gpu_packages <- paste0("cudatoolkit==", ver$cuda_version)
  }
  
  if (channel == "stable") {
    ver$channel <- "pytorch"
  } else if (channel == "nightly") {
    ver$channel <- "pytorch-nightly"
  } else {
    stop("not a valid channel")
  }
  
  ver
}



#' Install additional Python packages alongside PyTorch
#'
#' This function is deprecated. Use the `extra_packages` argument in function
#' `install_pytorch()` to install additional packages.
#'
#' @param packages Python packages to install
#' @param conda Path to conda executable (or "auto" to find conda using the PATH
#'   and other conventional install locations). Only used when PyTorch is
#'   installed within a conda environment.
#'
#' @keywords internal
#'
install_torch_extras <- function(packages, conda = "auto") {
  message("Extra packages not installed (this function is deprecated). \n",
          "Use the extra_packages argument to install_pytorch() to ",
          "install additional packages.")
}
