raw.fluorescence.check<-function(fcs.file.paths){
  keywords<-get.keywords(fcs.file.paths)
  if(any(keywords$parameters$TYPE %in% "Unmixed_Fluorescence")){
    stop(
      paste(
        "'fcs.files.paths' contains file(s) with 'Unmixed_Fluorescence' parameters.",
        "'fcs.files.paths' should contain only file paths to raw/reference .fcs files.",
        sep = "\n"
      )
    )
  }
}
##peak/cut function
peak.max.cut<-function(j){
  d<-stats::density(j)
  ##highest density peak index and value
  peak.i<-which.max(d$y)
  peak.x<-d$x[peak.i]
  ##peak cuts to select bounds
  low.cut<-mean(d$y[1:(peak.i-1)])
  low.cut.x<-d$x[d$y>low.cut][1]
  ##mirror the low cut
  high.cut.x<-peak.x+(peak.x-low.cut.x)
  return(c(low.cut.x,high.cut.x))
}
peaks.values<-function(j,quantile.trim=c(0.0001,0.9999),return.n.peaks=2){
  # d<-stats::density(j)
  q<-stats::quantile(j,probs = quantile.trim)
  d<-stats::density(j[j>q[1] & j<q[2]])
  peaks.i<-which(diff(sign(diff(d$y)))==-2)+1
  peaks.values<-d$x[peaks.i[order(d$y[peaks.i]) <= return.n.peaks]]
  return(peaks.values)
}
#' @title Get spectral signatures from raw reference controls
#' @description
#' Using Cytek Aurora raw reference controls -- .fcs files containing emission values for all available detectors -- either median values per detector or the raw values themselves can be returned.
#'
#' @param raw.fluor.control.paths Character string; path(s) usually returned from `list.files(...,full.names=T,pattern=".fcs")`. For Cytek Aurora experiments, the raw reference controls are stored in: `"./.../Raw/Reference Group"`.
#' @param return.type Character string; one of:
#' \itemize{
#'   \item `'medians'` -- return median expression values per-detector/per-control.
#'   \item `'raw'` -- return raw expression values per-detector/per-control.
#' }
#'
#' @returns a `spectr` object.
#' @export
#'
#' @examples
#'
#' extdata <- system.file("extdata",package="spectr")
#' raw.fluor.control.paths <- list.files(extdata,full.names=TRUE,pattern="Beads.*.fcs")
#'
#' spectral.signatures<-get.spectral.signature(raw.fluor.control.paths)
#' spectral.signatures$data[]
#'
#'
get.spectral.signature<-function(
    raw.fluor.control.paths,
    return.type=c('medians','raw')
){
  raw.fluorescence.check(raw.fluor.control.paths)
  dts<-read.fcs.to.dt(
    fcs.file.paths = raw.fluor.control.paths,
    colnames.syntatically.valid = TRUE,
    natural.order = FALSE,
    colnames.type = "N",
    cofactor = 5000,
    sample.id = "TUBENAME",
    keywords.to.factor = NULL
  )
  ##scatter variables; fluor variables
  # detectors.scatter<-dts$parameters[grep("scatter",TYPE,ignore.case = T),alias]
  detectors.fluors<-dts$parameters[TYPE=="Raw_Fluorescence",alias]
  ##density distributions to find singlet bead population;
  ##assumes singlet beads are represented by the highest density peak;
  ##cut FSC and SSC to isolate singlet beads using 'peak.max.cut' function
  dts$data[
    ,
    j = scatter.select :=
      data.table::`%between%`(FSC_A,peak.max.cut(FSC_A)) &
      data.table::`%between%`(SSC_A,peak.max.cut(SSC_A)),
    by = sample.id
  ]
  # ##test plot
  # ggplot2::ggplot(dts$data,ggplot2::aes(FSC_A,SSC_A)) +
  #   ggplot2::geom_hex(bins= 200) +
  #   ggplot2::facet_wrap(~sample.id) +
  #   ggplot2::geom_vline(
  #     data = dts$data[(scatter.select),.(lims=range(FSC_A)),by=sample.id],
  #     mapping = ggplot2::aes(xintercept = lims)
  #   ) +
  #   ggplot2::geom_hline(
  #     data = dts$data[(scatter.select),.(lims=range(SSC_A)),by=sample.id],
  #     mapping = ggplot2::aes(yintercept = lims)
  #   )
  # ##test counts
  # dts$data[,.N,keyby=.(scatter.select,sample.id)]
  ##mean to find the peak detector for each fluor and unstained
  detectors.peak<-dts$data[
    i = (scatter.select),
    j =.(detector.peak = names(which.max(sapply(.SD,mean)))),
    .SDcols = detectors.fluors,
    by = .(sample.id)
  ]
  ##factor based on detector order; laser line --> emission -->
  detectors.peak[,detector.peak:=factor(detector.peak,levels=detectors.fluors)]
  data.table::setorder(detectors.peak,detector.peak)
  detectors.peak[,sample.id := factor(sample.id,levels=as.character(sample.id))]
  dts$data[,sample.id := factor(sample.id,levels=levels(detectors.peak$sample.id))]
  ##test data
  dt.subset<-apply(detectors.peak,1,function(i){
    dts$data[
      i = (scatter.select) & sample.id==i[['sample.id']],
      j = .SD,
      .SDcols = i[['detector.peak']],
      by = sample.id]
  })
  dt.subset<-data.table::rbindlist(
    lapply(dt.subset,function(dt){
      dt[,detector:=dt[,names(.SD),.SDcols= is.numeric]]
      data.table::setnames(dt,dt[,unique(detector)],'value')
    })
  )
  dt.subset[,detector:=factor(detector,levels = detectors.fluors)]
  # ##test plot
  # ggplot2::ggplot(dt.subset,ggplot2::aes(value)) +
  #   ggplot2::geom_density() +
  #   ggplot2::facet_wrap(~sample.id,scales = "free")
  ##value of max peak (density) in peak detector using 'peak.max.x';
  ##exclude unstained
  for(i in detectors.peak[,seq(.N)]){
    if(!grepl("Unstained",detectors.peak[i,sample.id])){
      data.table::set(
        x = detectors.peak,
        i = i,
        j = paste0("detector.peak.",c('negative','positive')),
        value = as.list(dts$data[
          i = (scatter.select) & sample.id == detectors.peak[i,sample.id],
          j = sapply(.SD,peaks.values),
          .SDcols = detectors.peak[i,as.character(detector.peak)]
        ])
      )
    }
  }
  # ##test plot
  # ggplot2::ggplot(dt.subset,ggplot2::aes(value)) +
  #   ggplot2::geom_density() +
  #   ggplot2::geom_vline(
  #     data = data.table::melt(detectors.peak,id.vars=c('sample.id','detector.peak')),
  #     mapping = ggplot2::aes(xintercept = value)
  #   ) +
  #   ggplot2::facet_wrap(~sample.id,scales = "free")
  ##initialize 'events.select'; logical; FALSE
  ##set "Unstained" samples to TRUE
  ##bisect the max peak to select the top expressing events; fluors only
  dts$data[,events.select:=FALSE]
  dts$data[(scatter.select) & grepl("Unstained",sample.id),events.select := TRUE]
  for(i in detectors.peak[,seq(.N)]){
    data.table::set(
      x = dts$data,
      i = dts$data[
        ,
        j = .I[
          (scatter.select) &
            sample.id == detectors.peak[i,sample.id] &
            .SD>detectors.peak[i,detector.peak.positive]
        ],
        .SDcols = detectors.peak[i,as.character(detector.peak)]
      ],
      j = 'events.select',
      value = TRUE
    )
  }
  # ##test counts
  # dts$data[,.N,keyby=.(events.select,sample.id)]
  ##
  detectors.peak[,tissue.type := tolower(stringr::str_extract(sample.id,'Beads|Cells'))]
  detectors.peak[,marker.fluor:=sub(" \\(.*$","",sample.id)]
  detectors.peak[,marker:=sub(" .*","",marker.fluor)]
  detectors.peak[
    !marker %in% 'Unstained',
    fluorophore:=trimws(sub(paste0(marker,collapse = "|"),"",marker.fluor))
  ]
  detectors.peak[,c(paste0("detector.peak.",c('negative','positive')),'marker.fluor') := NULL]
  ##
  return.type<-match.arg(return.type)
  ##
  if(return.type=="raw"){
    spectr.transform.inverse(dts)
    dts$detectors.peak<-detectors.peak
    return(dts)
  }else if(return.type=='medians'){
    ##calculate medians/spectral signature
    detectors.medians<-dts$data[
      (events.select),
      lapply(.SD,stats::median),
      .SDcols = detectors.fluors,
      by=sample.id
    ]
    ##
    dts.medians<-c(
      list(data=detectors.medians[detectors.peak,on='sample.id']),
      dts[grep("data",attributes(dts)$names,invert = TRUE)],
      list(detectors.peak=detectors.peak)
    )
    spectr.transform.inverse(dts.medians)
    return(dts.medians)
  }
}
#' @title Plot a spectral signature
#' @description
#' Using Cytek Aurora raw reference controls -- .fcs files containing emission values for all available detectors -- median values per detector/per-control are calculated and then used to plot a 'trace' line representing a spectral signature.
#'
#' @inheritParams get.spectral.signature
#' @param subtract.universal Logical; default `TRUE`. Per-detector background fluorescence will be subtracted using a universal negative control, with respect to control-type (beads or cells).
#'
#' @returns A list containing \link[ggplot2]{ggplot2} objects.
#' @export
#'
#' @examples
#'
#' extdata <- system.file("extdata",package="spectr")
#' raw.fluor.control.paths <- list.files(extdata,full.names=TRUE,pattern="Beads.*.fcs")
#'
#' p<-spectral.signature.plot(raw.fluor.control.paths)
#' p[]
#'
spectral.signature.plot<-function(raw.fluor.control.paths,subtract.universal=TRUE){
  ##
  dts.medians<-get.spectral.signature(raw.fluor.control.paths,return.type = "medians")
  ##
  detectors.fluors<-dts.medians$parameters[TYPE=="Raw_Fluorescence",alias]
  ##
  if(subtract.universal){
    ##unstained medians; by 'tissue.type'
    medians.unstained.beads<-dts.medians$data[
      tissue.type == 'beads' & marker == "Unstained",
      .SD,
      .SDcols = detectors.fluors
    ]
    ##subtract unstained medians from fluors; by 'tissue.type'
    dts.medians$data[
      tissue.type == 'beads' & marker != "Unstained"
      ,
      (detectors.fluors):=unlist(.SD)-(medians.unstained.beads),
      .SDcols = detectors.fluors,
      by = sample.id
    ]
  }
  ##medians normalized to max median; by 'sample.id'
  dts.medians.normalized<-data.table::copy(dts.medians$data)[
    ,
    (detectors.fluors):=.SD/max(.SD),
    .SDcols = detectors.fluors,
    by = sample.id
  ]
  ##
  dt.melt<-data.table::melt(
    dts.medians.normalized,
    measure.vars=detectors.fluors
  )
  ##plot
  plot.list<-lapply(split(dt.melt,by='sample.id'),function(dt){
    ggplot2::ggplot(
      data=dt,
      mapping = ggplot2::aes(variable, value)
    ) +
      ggplot2::geom_line(linewidth=0.5,group=1) +
      ggplot2::geom_point() +
      ggplot2::geom_vline(
        xintercept = dt[,.I[variable==detector.peak]],
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
          droplevels(dt)[,levels(sample.id)]
        ),
        subtitle = paste(
          paste("Fluorophore:",dt[,unique(fluorophore)]),
          paste("Marker:", dt[,unique(marker)]),
          paste("Tissue Type:",dt[,unique(tissue.type)]),
          paste("Project:",unique(dts.medians$parameters.non$`$PROJ`)),
          sep = "\n"
        ),
        x = "Detector",
        y = "Normalized Median Expression",
        caption = paste("Peak Detector:", dt[,unique(detector.peak)])
      )
  })
  ##
  return(plot.list)
}

