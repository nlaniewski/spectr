titration.plot<-function(fcs.file.paths,detector.peak=NULL){
  dts<-read.fcs.to.dt(
    fcs.file.paths,
    colnames.syntatically.valid = TRUE,
    natural.order = FALSE,
    colnames.type = 'N',
    cofactor = 5000,
    sample.id = "TUBENAME"
  )
  ##
  dts$data[,c('test.volume','fluorophore','marker','clone'):=data.table::tstrsplit(sample.id,"_",type.convert = as.factor)]
  dts$data[,test.volume:=factor(as.numeric(sub("p",".",test.volume)))]
  ##
  scatter.cut.i<-dts$data[
    ,
    lapply(.SD,peak.max.cut),
    .SDcols = grep("[FS]SC",names(dts$data),value = T)
  ]

  scatter.cut.i<-Reduce(intersect,sapply(names(scatter.cut.i),function(j){
    dts$data[,.I[data.table::`%between%`(.SD,scatter.cut.i[[j]])],.SDcols = j]
  }))

  dts$data[,scatter.select := FALSE]
  dts$data[scatter.cut.i,scatter.select := TRUE]

  cols.fluors<-dts$parameters[TYPE=="Raw_Fluorescence",alias]

  if(is.null(detector.peak)){
    detector.peak<-names(which.max(table(dts$data[
      scatter.select==TRUE,
      names(which.max(unlist(lapply(.SD,stats::var)))),
      .SDcols = cols.fluors,
      by = sample.id][[2]])))
  }

  n<-dts$data[scatter.select==TRUE,.N,by = sample.id]
  .id.rand<-unlist(sapply(seq(nrow(n)),function(i){stats::rnorm(n=n[i,N],mean=10000*i,sd=1000)}))
  dts$data[scatter.select==TRUE,id.rand:=.id.rand]

  ##labels for plot
  fluor <- dts$data[,levels(fluorophore)]
  marker <- dts$data[,levels(marker)]
  clone <- dts$data[,levels(clone)]

  ggplot2::ggplot(
    dts$data[scatter.select==TRUE],
    ggplot2::aes(id.rand,!!as.name(detector.peak))
  ) +
    ggplot2::geom_hex(bins=300,show.legend = F) +
    viridis::scale_fill_viridis(option="plasma",limits = c(5,50),oob=scales::squish) +
    ggplot2::scale_x_continuous(
      breaks=seq(10000,n[,.N*10000],by=10000),
      labels=c(dts$data[,levels(test.volume)])
    ) +
    ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste("Titration:",fluor, marker, clone),
      # subtitle = (),
      x = "Test Volume (\u00B5L)",
      y = paste(
        detector.peak,
        fluor,
        marker,
        clone,
        sep = "::"
      )
    )
}
