#' @noRd
#' @title Get a list of all keyword-value pairs from .fcs file TEXT sections.
#' @param fcs.file.paths Character string; path(s) returned from `list.files(...,full.names=T,pattern=".fcs")`.
#' @return A list of length two; `[['parameters']]` and `[['parameters.non']]`.
get.keywords.text<-function(fcs.file.paths){
  ##read the .fcs text section; list (per-file); named character vector (all keywords)
  fcs.text.list<-flowCore::read.FCSheader(fcs.file.paths)
  ##keywords; list (per-file); named character vectors
  ##parameters ("$P|P"); non-parameters
  keywords<-lapply(fcs.text.list,function(fcs.text){
    sapply(c('parameters','parameters.non'),function(keyword.type){
      trimws(
        fcs.text[grep(
          pattern = "^\\$P\\d+|^P\\d+",
          x = names(fcs.text),
          invert = ifelse(keyword.type=='parameters.non',TRUE,FALSE))
        ]
      )
    },simplify = F)
  })
  ##
  return(keywords)
}
#' @noRd
#' @title Get a data.table of keyword parameters ("$P|P")
#' @param keywords.text The return of `get.keywords.text(...)`.
#' @param add.PROJ.identifier Logical; default `TRUE`. An additional column ('$PROJ') is added to the returned data.table.
get.dt.parameters<-function(keywords.text,add.PROJ.identifier=TRUE){
  #using the return of 'get.keywords.text(...)'
  lapply(keywords.text,function(kw){
    ##number of parameters
    par.n<-as.numeric(kw[['parameters.non']][['$PAR']])
    ##parameters vector
    pars<-kw[['parameters']]
    ##split parameters vector into 'types'; list
    par.types<-split(x=pars,f=factor(sub("^\\$P\\d+|^P\\d+","",names(pars))))
    ##create a data.table
    dt.parameters<-data.table::data.table(par=paste0('$P',seq(par.n)))
    ##fill the data.table by 'par.type'
    for(j in names(par.types)){
      par.vec<-par.types[[j]]
      par.i<-as.integer(gsub("\\D+","",names(par.vec)))
      data.table::set(dt.parameters,i=par.i,j=j,value = par.vec)
    }
    ##add '$PROJ' identifier
    if(add.PROJ.identifier){
      dt.parameters[,PROJ:=as.factor(kw[['parameters.non']][['$PROJ']])]
    }
    ##return the data.table
    dt.parameters[]
  })
}
#' @noRd
#' @title Get a data.table of non-parameter keywords.
#' @param keywords.text The return of `get.keywords.text(...)`.
#' @param drop.primary Logical; default `TRUE`. Drops required keyword-value pairs from the returned data.table.
#' @param drop.spill Logical; default `TRUE`. Drops the lengthy spillover string from the returned data.table.
get.dt.parameters.non<-function(keywords.text,drop.primary=TRUE,drop.spill=TRUE){
  #using the return of 'get.keywords.text(...)'
  lapply(keywords.text,function(kw){
    ##'parameters.non' vector to data.table
    dt.parameters.non<-data.table::as.data.table(as.list(kw[['parameters.non']]))
    ##drop required keywords
    if(drop.primary){
      dt.parameters.non[,fcs.text.primary.required.keywords():=NULL]
    }
    if(drop.spill){
      spill.name<-grep("spill",names(dt.parameters.non),ignore.case = T,value = T)
      dt.parameters.non[,(spill.name):=NULL]
    }
    ##return the data.table
    dt.parameters.non[]
  })
}
get.spill<-function(keywords.text,add.PROJ.identifier=TRUE){
  #using the return of 'get.keywords.text(...)'
  lapply(keywords.text,function(kw){
    spill.name<-grep('spill',names(kw[['parameters.non']]),ignore.case = T,value = T)
    spillover.string<-kw[['parameters.non']][[spill.name]]
    ##from flowCore:::txt2spillmatrix
    spill.split <- strsplit(spillover.string, ",")[[1]]
    N.cols <- as.numeric(spill.split[1])
    if(!is.na(N.cols)&&N.cols>0){
      col.names <- spill.split[2:(N.cols + 1)]
      vals<-as.numeric(spill.split[(N.cols+2):length(spill.split)])
      spill.mat<-matrix(
        data = vals,
        ncol = N.cols,
        byrow = TRUE,
        dimnames = list(NULL,col.names)
      )
      if(add.PROJ.identifier){
        attr(spill.mat,'PROJ')<-kw[['parameters.non']][['$PROJ']]
        return(spill.mat)
      }
      return(spill.mat)
    }else{
      stop("'spillover.string' cannot be parsed.")
    }
  })
}
#' @title Get keywords from .fcs file(s)
#' @description
#' Parameter ("$P|P") keywords in the TEXT section of a .fcs file specify standard, cytometer-specific parameter names/values; Non-parameter keywords in the TEXT section of a .fcs file specify standard/optional cytometer/sample-specific meta-data. The "$SPILLOVER" keyword in the TEXT section of a .fcs file specifies a spillover matrix (for fluorescence compensation). This function allows parsing of these keywords without the overhead of reading the DATA section of the .fcs file.
#'
#' @param fcs.file.paths Character string; path(s) usually returned from `list.files(...,full.names=T,pattern=".fcs")`.
#'
#' @returns a list containing three elements:
#' \itemize{
#'    \item `[['parameters']]` -- parameters, unique by 'PROJ' (from '$PROJ': name(s) of the experiment project(s)); row-bound `data.table`.
#'    \item `[['parameters.non']]` -- non-parameter keywords; row-bound `data.table`.
#'    \item `[['spill']]` -- a list of spillover matrices; unique by 'PROJ' (from '$PROJ': name(s) of the experiment project(s)).
#' }
#'
#' @export
#'
#' @examples
#'
#' extdata <- system.file("extdata",package="spectr")
#' fcs.file.paths <- list.files(extdata,full.names=TRUE,pattern="BLOCK.*.fcs")
#'
#' keywords<-get.keywords(fcs.file.paths)
#'
#' keywords$parameters #parameters data.table
#' keywords$parameters.non #non-parameters data.table
#' keywords$spill #spill/spillover
#'
#' #project-specific parameters
#' keywords$parameters[PROJ %in% levels(PROJ)[1]]
#' #project-specific parameters; selecting two columns; data.table syntax
#' keywords$parameters[PROJ %in% levels(PROJ)[1]][,.(N,S)]
#'
#' #project-specific keywords
#' keywords$parameters.non[`$PROJ` %in% unique(`$PROJ`)[1]]
#' #project-specific keywords; selecting two columns; data.table syntax
#' keywords$parameters.non[`$PROJ` %in% unique(`$PROJ`)[1]][,.(TUBENAME,`$DATE`)]
#'
#' #create new keywords
#' keywords$parameters.non[
#' i=grep("001_CYTOKINE",TUBENAME),
#' j=c('study.id','batch.id','panel.id','block.id','block.aliquot')
#' :=data.table::tstrsplit(TUBENAME,"_",type.convert=as.factor,keep=c(1:3,5:6))]
#'
#' #create new keywords
#' keywords$parameters.non[
#' i=grep("AIM|002_CYTOKINE",TUBENAME),
#' j=c('study.id','batch.id','panel.id','block.id','block.aliquot')
#' :=data.table::tstrsplit(TUBENAME,"_",type.convert=as.factor)]
#'
#' #construct a new 'sample.id' from existing keywords;
#' #use data.table's assignment operator ':='
#' keywords$parameters.non[,
#' sample.id := paste(study.id,batch.id,panel.id,block.id,block.aliquot,sep="_")]
#' keywords$parameters.non[['sample.id']]
#'
get.keywords<-function(fcs.file.paths){
  ##keywords for each file; list
  keywords.text<-get.keywords.text(fcs.file.paths)
  ##parameters for each file; list
  dts.parameters<-get.dt.parameters(keywords.text)
  ##row-bind
  dt.parameters<-data.table::rbindlist(dts.parameters,fill=TRUE)
  ##unique dt.parameters; by 'PROJ'
  ##probably an easier way, but below code chunk works to find unique 'Time' by max 'R'
  dt.parameters<-data.table::rbindlist(
    lapply(split(dt.parameters,by='PROJ'),function(dt.proj){
      dt<-unique(dt.proj)
      non.unique.par<-dt[,which(.N>1),by=par][,par]
      NR<-unique(dt[par %in% non.unique.par,if(length(unique(N))==1){.(N,R=max(R))}])
      dt[-dt[,.I[N %in% NR[['N']]&R!=NR[['R']]]]][order(as.integer(gsub("\\D+","",par)))]
    })
  )
  ##non-parameter keywords for each file; list
  dts.parameters.non<-get.dt.parameters.non(keywords.text)
  ##row-bind
  dt.parameters.non<-data.table::rbindlist(dts.parameters.non,fill=TRUE)
  ##spill for each file; list
  ##unique 'spill.mats'; by 'PROJ'
  spill.mats<-unique(get.spill(keywords.text))
  ##return
  list(
    parameters=dt.parameters,
    parameters.non=dt.parameters.non,
    spill=spill.mats
  )
}
