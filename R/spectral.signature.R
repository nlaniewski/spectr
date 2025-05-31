#' @title Plot a normalized spectral signature
#'
#' @description
#' Plot a normalized spectral signature for spectral flow cytometry reference controls.
#'
#' @param raw.fluor.control.path Character string; file path to a spectral flow cytometry reference (fluor) control.
#' @param raw.unstained.control.path Character string; file path to a spectral flow cytometry reference (unstained, universal) control to be used for signal (auto fluorescence) subtraction. If `NULL`, an internal negative/unstained population will used.
#' @param transform.type Character string; default and only current choice is 'asinh' for arc-sine transformation of fluorescence detectors.
#' @param cofactor Numeric; default 5000. Defines the cofactor value used during arc-sine transformation of fluorescence detectors.
#' @param top.n Numeric; default 200. Defines the top number of expressing events (peak detector) to use when generating the normalized median values.
#'
#' @return A \link[ggplot2]{ggplot} object.
#' @export
#' @references https://github.com/hally166/flowSpectrum/
#'
#' @examples
#' raw.fluor.control.paths<-list.files(
#' system.file("extdata",package="spectr"),
#' full.names = TRUE,
#' pattern = "RB744|PerCPFire806"
#' )
#'
#' raw.unstained.control.path<-list.files(
#' system.file("extdata",package="spectr"),
#' full.names = TRUE,
#' pattern = "Unstained"
#' )
#'
#' ##spectral signatures using universal negative/unstained for signal subtraction
#' p.universal<-lapply(
#' raw.fluor.control.paths,
#' spectral.signature,
#' raw.unstained.control.path
#' )
#'
#' ##spectral signatures using internal negative/unstained for signal subtraction
#' p.internal<-lapply(
#' raw.fluor.control.paths,
#' spectral.signature,
#' raw.unstained.control.path=NULL
#' )
#'
#' ##CD3 RB744
#' p.universal[[1]] +
#' ggplot2::labs(caption="Universal negative/unstained used for signal subtraction.")
#' p.internal[[1]] +
#' ggplot2::labs(caption="Internal negative/unstained used for signal subtraction.")
#'
#' ##CD4 PerCPFire806
#' p.universal[[2]] +
#' ggplot2::labs(caption="Universal negative/unstained used for signal subtraction.")
#' p.internal[[2]] +
#' ggplot2::labs(caption="Internal negative/unstained used for signal subtraction.")
#'
#'
spectral.signature<-function(
    raw.fluor.control.path,
    raw.unstained.control.path,
    transform.type='asinh',
    cofactor=5000,
    top.n=200){
  ##to silence NSE R CMD check notes
  FSC_A<-SSC_A<-detector<-normalized.expression<-NULL

  ##https://github.com/hally166/flowSpectrum/
  ##normalized spectral signatures

  ##get .fcs text/keywords
  fcs.kw<-as.list(flowCore::read.FCSheader(raw.fluor.control.path)[[1]])
  fcs.kw.meta<-fcs.kw[grep('CYT|CYTSN|CREATOR|PROJ|TUBENAME',names(fcs.kw))]
  if(!any(unlist(unique(fcs.kw[grep("\\$P[0-9]+TYPE",names(fcs.kw))])) %in% "Raw_Fluorescence")){
    stop("'Raw_Fluorescence' ($P#TYPE) channels not detected; is this a reference control file (raw not unmixed)?")
  }
  ##list file paths
  dts<-list(
    'fluor'=raw.fluor.control.path,
    'unstained'=if(!is.null(raw.unstained.control.path)){
      raw.unstained.control.path
    }
  )
  ##if raw.unstained.control.path is not-defined, drop
  dts<-dts[!sapply(dts,is.null)]
  ##read .fcs data; convert to data.table
  dts<-lapply(dts,function(fcs.path){
    data.table::as.data.table(
      flowCore::read.FCS(
        fcs.path,
        transformation = F,truncate_max_range = F
      )@exprs
    )
  })
  ##change column names; syntactically valid names
  old.names<-data.table::copy(names(dts$fluor))
  if(!all(old.names==names(dts$unstained))){
    stop("Column names do not match between unstained control and raw control.")
  }
  new.names<-sub("SSC-B","SSCB",old.names)
  new.names[grep("[FS]SC",new.names)]<-sub("-","_",new.names[grep("[FS]SC",new.names)])
  new.names<-sub("-A$","",new.names)
  ##update dts by reference
  invisible(lapply(dts,function(dt){data.table::setnames(dt,old.names,new.names)}))
  ##scatter variables; fluor variables
  detectors.scatter<-grep("[FS]SC",names(dts$fluor),value = T)
  detectors.fluors<-grep("Time|[FS]SC",names(dts$fluor),value = T,invert = T)
  ##transform the fluors; updates by reference
  invisible(lapply(dts,function(dt){
    if(transform.type=='asinh'){
      for(j in detectors.fluors){
        data.table::set(dt,j=j,value=asinh(dt[[j]]/cofactor))
      }
    }
  }))
  ##density distributions to agnostically find singlet bead population
  ##assumes singlet beads are represented by the max peak
  ##peak/cut function
  peak.cut<-function(j){
    d<-stats::density(j)
    ##peak index and value
    peak.i<-which.max(d$y)
    peak.x<-d$x[peak.i]
    ##peak cuts to select bounds
    low.cut<-mean(d$y[1:(peak.i-1)])
    low.cut.x<-d$x[d$y>low.cut][1]
    ##mirror the low cut
    high.cut.x<-peak.x+(peak.x-low.cut.x)
    return(c(low.cut.x,high.cut.x))
  }
  ##cut FSC and SSC to isolate singlet beads
  dts<-lapply(dts,function(dt){
    scatter.cuts<-dt[,lapply(.SD,peak.cut),.SDcols = c('FSC_A','SSC_A')]
    dt<-dt[data.table::`%between%`(FSC_A,scatter.cuts$FSC_A)&data.table::`%between%`(SSC_A,scatter.cuts$SSC_A)]
    return(dt)
  })
  ##peak detector for fluor control
  detector.means<-dts$fluor[,lapply(.SD,function(x){
    mean(sort(x,decreasing = T)[1:top.n])
  }),.SDcols = detectors.fluors]
  detector.max<-names(which.max(detector.means))
  ##select top number of events
  detector.index<-dts$fluor[,.I[get(detector.max)>=sort(get(detector.max),decreasing = T)[top.n]]]
  ##per-detector medians; top.n events; fluor
  dt.medians<-dts$fluor[detector.index,lapply(.SD,stats::median),.SDcols = detectors.fluors]
  ##per-detector medians; negative/unstained; to subtract (auto fluorescence)
  if(!is.null(raw.unstained.control.path)){
    ##per-detector medians; unstained; universal
    dt.medians.unstained<-dts$unstained[,lapply(.SD,stats::median),.SDcols = detectors.fluors]
  }else{
    ##per-detector medians; negative/unstained; internal
    ##200 events centered over negative/unstained peak
    d<-dts$fluor[,stats::density(get(detector.max))]

    d.derivative.1<-diff(d$y)
    d.derivative.2<-diff(sign(d.derivative.1))
    peaks<-which(d.derivative.2==(-2))+1
    peak1<-d$x[peaks[1]]
    bound1<-dts$fluor[,sort(get(detector.max)[get(detector.max)<peak1],decreasing = T)[(top.n/2)]]
    bound2<-dts$fluor[,sort(get(detector.max)[get(detector.max)>peak1],decreasing = F)[(top.n/2)]]

    # plot(d)
    # abline(v=peak1,col="blue"))

    dt.medians.unstained<-dts$fluor[data.table::between(get(detector.max),bound1,bound2),
                                    lapply(.SD,stats::median),
                                    .SDcols = detectors.fluors]
  }
  ##subtract unstained medians from fluor
  dt.medians<-(dt.medians)-(dt.medians.unstained)
  ##medians normalized to max median
  dt.medians.normalized<-data.table::data.table(
    detector=names(dt.medians),normalized.expression=unlist(dt.medians/max(dt.medians))
  )
  ##factor the 64 detectors: UV -> V -> B -> YG -> R
  dt.medians.normalized[,detector := factor(detector,levels = detector)]
  ##plot
  ggplot2::ggplot(
    data=dt.medians.normalized,
    mapping = ggplot2::aes(detector, normalized.expression)
  ) +
    ggplot2::geom_line(linewidth=0.5,group=1) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(
      xintercept = dt.medians.normalized[,.I[detector==detector.max]],
      linetype = 'dashed'
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        color = "black",
        angle = 90,
        vjust = 0.5,
        hjust=1)
    ) +
    ggplot2::labs(
      title = paste0(
        "Spectral Signature: Normalized Median Expression",
        paste0(rep(" ",10),collapse = ""),
        fcs.kw.meta$TUBENAME
      ),
      # subtitle = paste(
      #   paste0(names(marker.fluor.type.batch),':'),
      #   gsub("\\(|\\)","",marker.fluor.type.batch),
      #   collapse = "\n"
      # ),
      subtitle = paste(
        paste0(names(fcs.kw.meta),':'),
        fcs.kw.meta,
        collapse = "\n"
      ),
      x = "Detector",
      y = "Normalized Median Expression",
      caption = paste("Peak Detector:", detector.max)
    )
}
