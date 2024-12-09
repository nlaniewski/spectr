.datatable.aware = TRUE
#for R CMD check; data.table vars
utils::globalVariables(
  c(
    '.I'
    ,'.SD'
    ,':='
    ,'FSC_A'
    ,'SSC_A'
    ,'detector'
    ,'normalized.expression'
  )
)
