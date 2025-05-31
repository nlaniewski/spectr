.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(strsplit(trimws(utils::readClipboard())," ")[[1]])
utils::globalVariables(
  c("$TOT", ".", "N", "N.alias", "PROJ", "R", "S", "S_N", "TYPE",
    "alias", "channels", "get.keywords.dt", "par", "spill.from.string")
)
