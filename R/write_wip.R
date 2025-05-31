# devtools::document();devtools::load_all()
# extdata <- system.file("extdata",package="spectr")
# fcs.file.paths <- list.files(extdata,full.names=TRUE,pattern="CYTOKINE.*BLOCK.*.fcs")
#
# dts<-read.fcs.to.dt(fcs.file.paths,keywords.to.factor = c('TUBENAME','$DATE'))
#
# cols.factor<-dts[,names(.SD),.SDcols = is.factor]
# if('TUBENAME' %in% cols.factor){split.by<-'TUBENAME'}
#
# lapply(split(dts,by=split.by),function(dt){
#   id<-dt[,levels(droplevels(j)),env = list(j=split.by)]
#   attr(dt,'keywords')[i==id,env = list(i=split.by)]
# })
#
#
#
# View(SOMnambulate:::parms.list.from.adf)
# View(SOMnambulate:::parms.adf.dt)
#
# fcs<-flowCore::flowFrame(
#   exprs = as.matrix(dt[TUBENAME==levels(TUBENAME)[1],-'TUBENAME'])
# )
