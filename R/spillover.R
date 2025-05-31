# spill.from.string<-function(spillover.string,channel_alias=NULL){
#   # spillover.string<-unlist(flowCore::read.FCSheader(fcs.file.path,keyword = '$SPILLOVER'))[[1]]
#   ##from flowCore:::txt2spillmatrix
#   spill.split <- strsplit(spillover.string, ",")[[1]]
#   N.cols <- as.numeric(spill.split[1])
#   if(!is.na(N.cols)&&N.cols>0){
#     col.names <- spill.split[2:(N.cols + 1)]
#     vals<-as.numeric(spill.split[(N.cols+2):length(spill.split)])
#     spill.mat<-matrix(
#       data = vals,
#       ncol = N.cols,
#       byrow = TRUE,
#       dimnames = list(NULL,col.names)
#     )
#   }else{
#     stop("'spillover.string' cannot be parsed.")
#   }
#   if(!is.null(channel_alias)){
#     if(!all(colnames(spill.mat) %in% channel_alias$channels)){
#       stop("Column name mismatch between 'spillover.string' and 'channel_alias$channels'.")
#     }else{
#       colnames(spill.mat)<-channel_alias$alias[match(colnames(spill.mat),channel_alias$channels)]
#     }
#   }
#   ##
#   return(spill.mat)
# }
# spill.update.value<-function(spill,i,j,value){
#   i.grep<-grep(i,colnames(spill))
#   if(length(i.grep)!=1){
#     stop(sprintf("grep for %s returned %s values; modify the 'pattern' (i) argument",
#                  i,length(i.grep)))
#   }
#   j.grep<-grep(j,colnames(spill))
#   if(length(j.grep)!=1){
#     stop(sprintf("grep for %s returned %s values; modify the 'pattern' (i) argument",
#                  j,length(j.grep)))
#   }
#   spill[i.grep,j.grep]<-value
#   return(spill)
# }
