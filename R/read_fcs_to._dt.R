concatenate.test<-function(fcs.file.paths){
  keywords<-get.keywords(fcs.file.paths)
  dt.parameters.split<-split(keywords[['parameters']],by='PROJ')
  dt.NS<-unique(
    lapply(dt.parameters.split,function(dt.parameters){dt.parameters[,.(N,S)]})
  )
  if(length(dt.NS)!=1){
    discrepancy<-sapply(c('N','S'),function(j){
      j.intersect<-Reduce(intersect,lapply(dt.NS,'[[',j))
      j.union<-Reduce(union,lapply(dt.NS,'[[',j))
      j.discrepancy<-j.union[!j.union %in% j.intersect]
    },simplify = F)
    discrepancy<-discrepancy[sapply(discrepancy,length)>1]
    ##
    stop(
      paste(
        sprintf("'$P%s' discrepancy; %s conflicts found: %s",
                names(discrepancy),
                sapply(discrepancy,length),
                sapply(discrepancy,paste,collapse=", ")
        ),
        "These files cannot be concatenated; are they all from the same project/staining panel?",
        paste(
          "To resolve (if you are certain they are from the same panel; maybe there was a typo?):",
          "1) set the argument 'return.list' to TRUE;",
          "2) name-correct list items (data.tables);",
          "3) 'data.table::rbindlist(...)'",
          sep = "\n"
        ),
        sep = "\n"
      )
    )
  }else{
    dt.NS<-dt.NS[[1]]
  }
  if(dt.NS[!is.na(S),data.table::uniqueN(S)!=.N]){
    stop("Non-unique '$PS' names")
  }else{
    return(keywords)
  }
}
prepare.channel.alias<-function(
    dt.parameters,
    colnames.syntatically.valid=TRUE,
    natural.order=TRUE
){
  ##to silence NSE R CMD check notes

  ##
  channel_alias<-unique(dt.parameters[,.(N,S)])
  ##
  if(colnames.syntatically.valid){
    channel_alias[!grepl("[FS]SC|Time",N),N.alias:=gsub(" |-","",sub("-A$","",N))]
    channel_alias[grepl("[FS]SC",N),N.alias:=sub("-","_",sub("SSC-B","SSCB",N))]
  }else{
    channel_alias[,N.alias:=N]
  }
  ##
  channel_alias[!is.na(S),S_N:=paste(S,N.alias,sep = "_")]
  channel_alias[is.na(S_N),S_N:=N.alias]
  channel_alias[is.na(S_N),S_N:=N]
  ##
  if(natural.order){
    channel_alias<-channel_alias[stringr::str_order(S_N,numeric = T)]
  }
  ##
  channel_alias[]
}
#' @title Read .fcs files and concatenate (row-bind) into a \link[data.table]{data.table}.
#'
#' @param fcs.file.paths Character string; path(s) usually returned from `list.files(...,full.names=T,pattern=".fcs")`.
#' @param colnames.syntatically.valid Logical; default `TRUE`. Syntactically valid names (letters, numbers and the underline character) are enforced by removing all dashes and spaces from column/parameter names.
#' @param natural.order Logical; default `TRUE`. Naturally order column names.
#' @param colnames.type Character string; one of:
#' \itemize{
#'   \item `"S_N"` -- parameter columns are named by combining $PS and $PN, separated by an underscore. For Cytek Aurora instruments using SpectroFlo software, the default parameter naming scheme is as follows: $PS = stain (marker); $PN = detector/fluorophore.
#'   \item `"S"` -- parameter columns are named by using only their respective $PS keyword value.
#' }
#' @param cofactor Numeric; default `5000`. Any/all parameters with a `$PnTYPE` of 'Unmixed_Fluorescence' will be transformed using \link{asinh} and the defined cofactor value (`asinh(x/cofactor)`).
#' @param sample.id Character string; the keyword label defined through `sample.id` (default `TUBENAME`) will be used to add respective keyword values as an identifier to the `data.table`.
#' @param keywords.to.factor Character string; if defined, keywords will be coerced to factor and appended to the returned data.table as factored columns.
#'
#' @returns a list
#' @export
#'
#' @examples
#'
#' extdata <- system.file("extdata",package="spectr")
#' fcs.file.paths <- list.files(extdata,full.names=TRUE,pattern="BLOCK.*.fcs")
#'
#' #trying to read mixed panels ("AIM","CYTOKINE") will cause an error due to $PS name conflicts
#' \dontrun{
#' read.fcs.to.dt(fcs.file.paths)
#' }
#'
#' #"CYTOKINE" panel
#' fcs.file.paths<-grep("CYTOKINE",fcs.file.paths,value=TRUE)
#'
#' dt<-read.fcs.to.dt(fcs.file.paths)
#'
#' names(dt)
#' sapply(c("data","parameters","parameters.non"),function(i){head(dt[[i]])})
#' head(dt$spill[[1]])
#'
#' #add additional (factored) keywords
#' dt<-read.fcs.to.dt(fcs.file.paths,keywords.to.factor=c('$PROJ','$DATE'))
#'
#' dt$data
#'
read.fcs.to.dt<-function(
    fcs.file.paths,
    colnames.syntatically.valid=TRUE,
    natural.order=TRUE,
    colnames.type=c("S","S_N"),
    cofactor=5000,
    sample.id='TUBENAME',
    keywords.to.factor=NULL
    # return.list=FALSE,
){
  ##to silence NSE R CMD check notes

  ##
  keywords<-concatenate.test(fcs.file.paths)
  ##
  channel_alias<-prepare.channel.alias(
    keywords$parameters,
    colnames.syntatically.valid,
    natural.order
  )
  ##
  ctype<-match.arg(colnames.type)
  if(ctype=='S_N'){
    ca<-channel_alias[,.(channels=N,alias=S_N)]
  }else if(ctype=='S'){
    ca<-channel_alias[,.(channels=N,alias=S)]
    ca[is.na(alias),alias:=channel_alias[N %in% ca[is.na(alias),channels],S_N]]
  }
  ##
  keywords$parameters<-merge(keywords$parameters,ca[,.(N=channels,alias)],sort = F)
  ##
  dt<-data.table::rbindlist(
    lapply(fcs.file.paths,function(fcs.file.path){
      data.table::as.data.table(
        flowCore::read.FCS(
          filename = fcs.file.path,
          transformation = F,
          truncate_max_range = F,
          channel_alias = ca
        )@exprs
      )
    })
  )
  ##
  if(natural.order){
    data.table::setcolorder(dt, ca$alias)
  }
  ##
  if(!is.null(cofactor)){
    cols.transform<-keywords$parameters[TYPE=="Unmixed_Fluorescence",unique(alias)]
    for(j in cols.transform){
      data.table::set(dt,j=j,value = asinh(dt[[j]]/cofactor))
    }
  }
  ##
  data.table::set(
    dt,
    j='sample.id',
    value=as.factor(
      keywords$parameters.non[,rep(j,as.numeric(`$TOT`)),
        env = list(j = sample.id)])
  )
  ##
  if(!is.null(keywords.to.factor)){
    for(i in seq(keywords.to.factor)){
      data.table::set(
        dt,
        j = sub("\\$","",keywords.to.factor[i]),
        value=as.factor(
          keywords$parameters.non[,rep(j,as.numeric(`$TOT`)),
                                  env = list(j = keywords.to.factor[i])])
      )
    }
  }
  ##
  return(c(list(data=dt),keywords))
  ##
}
