raw.fluor.control.paths<-c(
  "./inst/extdata/CD3_RB744_Beads.fcs",
  "./inst/extdata/CD4_PerCPFire806_Beads.fcs",
  "./inst/extdata/Unstained_Beads.fcs"
)

fcs.file.paths<-c(raw.fluor.control.paths,"./inst/extdata/COVAIL_001_CYTOKINE_BCpool_BLOCK1_1.fcs")

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

peak.max.x<-function(j){
  d<-stats::density(j)
  peaks.i<-which(diff(sign(diff(d$y)))==-2)+1
  peak.x.max<-d$x[max(peaks.i)]
}

spectral.signature.medians<-function(
    raw.fluor.control.paths
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
  ##cut FSC and SSC to isolate singlet beads
  scatter.cuts<-dts$data[,lapply(.SD,peak.max.cut),.SDcols = c('FSC_A','SSC_A'),by=sample.id]
  # ggplot2::ggplot(
  #   data = dts$data[,.(FSC_A,SSC_A),by=sample.id],
  #   mapping = ggplot2::aes(FSC_A,SSC_A,fill=sample.id)
  # ) +
  #   ggplot2::geom_hex(bins=200) +
  #   ggplot2::geom_vline(xintercept = scatter.cuts$FSC_A) +
  #   ggplot2::geom_hline(yintercept = scatter.cuts$SSC_A)
  dts$data[,scatter.select:=FALSE]
  for(s in scatter.cuts[,levels(sample.id)]){
    data.table::set(
      dts$data,
      i = dts$data[
        ,
        .I[sample.id == s &
             data.table::`%between%`(FSC_A,scatter.cuts[sample.id == s,FSC_A]) &
             data.table::`%between%`(SSC_A,scatter.cuts[sample.id == s,SSC_A])]
      ],
      j = 'scatter.select',
      value = TRUE
    )
  }
  # dts$data[,.N,by=.(scatter.select,sample.id)]
  ##variance to find the peak detector for each fluor
  detectors.max<-dts$data[
    scatter.select==TRUE,
    .(detector=names(which.max(unlist(sapply(.SD,var))))),
    .SDcols = detectors.fluors,
    by=sample.id
  ]
  ##value of max peak (density) in peak detector
  detectors.max[,peak.max.value:=0]
  for(s in detectors.max[,levels(sample.id)]){
    data.table::set(
      x = detectors.max,
      i = detectors.max[,.I[sample.id == s]],
      j = 'peak.max.value',
      value = dts$data[
        sample.id == s & scatter.select == TRUE,
        peak.max.x(j),
        env = list(j = detectors.max[sample.id == s,detector])
      ]
    )
  }
  ##bisect the max peak to select the top expressing events
  dts$data[,events.select:=FALSE]
  for(s in detectors.max[,levels(sample.id)]){
    data.table::set(
      dts$data,
      i = dts$data[
        ,
        .I[sample.id == s & scatter.select == TRUE & d > val],
        env = list(
          d = detectors.max[sample.id == s, detector],
          val = detectors.max[sample.id == s, peak.max.value]
        )
      ],
      j = 'events.select',
      value = TRUE
    )
  }
  # dts$data[,.N,by=.(events.select,sample.id)]
  ##calculate medians/spectral signature
  detectors.medians<-dts$data[
    events.select==TRUE,
    lapply(.SD,stats::median),
    .SDcols = detectors.fluors,
    by=sample.id
  ]
  ##
  return(
    c(
      list(data=detectors.medians),
      dts[grep("data",attributes(dts)$names,invert = TRUE)],
      list(detectors.max=detectors.max[,-'peak.max.value'])
    )
  )
}

