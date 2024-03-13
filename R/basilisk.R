#' @title xxx
#' @description xxx
#' @importFrom basilisk BasiliskEnvironment
#' @export
#' @return An instance of the class `BasiliskEnvironment`
#' 
envPirat <- basilisk::BasiliskEnvironment("envPirat",
  pkgname="Pirat",
  packages=c("numpy==1.20.2",
    "matplotlib==3.7",
    "pytorch==1.10.0",
    "python==3.9.5"),
  channels = c("pytorch", "stable", "torch"),
  path="test_dummy")
