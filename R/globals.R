.datatable.aware = TRUE
##to silence NSE R CMD check notes; "no visible biding for global variable..."
# dput(unlist(strsplit(trimws(utils::readClipboard())," ")))
utils::globalVariables(
  c(
    c("$TOT", ".", "N", "N.alias", "PROJ", "R", "S", "S_N", "TYPE",
      "alias", "channels", "get.keywords.dt", "par", "spill.from.string"),
    c("FSC_A", "SSC_A", "cofactor", "detector", "detector.peak",
      "detector.peak.positive", "events.select", "fluorophore", "marker",
      "marker.fluor", "quantile", "sample.id", "scatter.select", "tissue.type",
      "value", "variable")
  )
)
