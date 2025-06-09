spectr.transform.inverse<-function(spectr.object){
  dt.transform<-spectr.object$parameters[grep("Raw|Unmixed_Fluorescence",TYPE),.(alias,cofactor)]
  for(j in dt.transform$alias){
    data.table::set(
      x = spectr.object$data,
      j = j,
      value = sinh(spectr.object$data[[j]])*dt.transform[alias==j,cofactor]
    )
  }
  spectr.object$parameters[,c('transform','cofactor') := NULL]
}
