fcs.synthetic.tmp.files<-function(
    n.files=9,
    events=10,
    parameters=26,
    error.type=c('none','extra keyword','$PN conflict','$PS conflict')
){
  if(!requireNamespace("flowCore")){
    stop("This function needs the flowCore package.")
  }
  if(n.files>9){
    stop("'n.files' value cannot exceed 9")
  }
  if(parameters>26){
    stop("'parameters' value cannot exceed 26")
  }
  for(i in seq(n.files)){
    ##synthetic data
    mat1<-matrix(
      data = {set.seed(i);stats::rnorm(events*parameters)},
      nrow = events,
      ncol = parameters,
      dimnames = list(NULL,paste0("detector ",LETTERS[seq(parameters)],"-A"))
    )
    ##synthetic data -- time and scatter
    mat2<-matrix(
      data = c(sort(sample.int(n=(events/5000/0.0001),size = events,replace = TRUE)),
               sample(0:4194304,events),
               sample(0:4194304,events)),
      nrow = events,
      ncol = 3,
      dimnames = list(NULL,c("Time","FSC-A","SSC-A"))
    )
    mat<-cbind(mat2,mat1)
    ##optional keyword; 9th keyword sample-specific
    animals<-c("aardvark","beaver","capybara","dog","elephant","frog","gecko","hamster","impala")
    keyword.string<-sprintf("The quick brown fox jumps over the lazy %s",animals[i])
    keyword.vec<-stats::setNames(
      strsplit(keyword.string," ")[[1]],
      nm=paste0("KEYWORD",1:9)
    )
    ##initialize keywords (list)
    keywords<-as.list(keyword.vec)
    ##additional keywords; standard ('$') and optional;
    ##following FSC 3.1 normative reference
    keywords[['$BTIM']]<-data.table::as.ITime("12:00:00.00")+(i*10)#hh:mm:ss[.cc]; data.table class ("ITime") gets stripped when writing .fcs file
    keywords[['$COM']]<-"This .fcs file was built using R and is for example/testing purposes."
    keywords[['$CYT']]<-"InSilicometer"
    keywords[['$CYTSN']]<-"1337"
    keywords[['$DATE']]<-toupper(format(Sys.Date(),"%d-%b-%Y"))
    keywords[['$ETIM']]<-keywords[['$BTIM']]+10#add 10 seconds
    # keywords[['FCSversion']]<-3.1#flowCore::write.FCS produces version 3.0
    keywords[['$FIL']]<-sprintf("sample%s.fcs",LETTERS[i])
    keywords[['$INST']]<-"URMC Schieble Lab"
    keywords[['$LAST_MODIFIED']]<-{options(digits.secs=2);toupper(format(Sys.time(), "%d-%b-%Y %H:%M:%OS"))}
    keywords[['$LAST_MODIFIER']]<-"Nathan Laniewski"
    keywords[['$ORIGINALITY']]<-"DataModified"
    keywords[['$PROJ']]<-"spectr_package"
    keywords[['sample.id']]<-sprintf("sample%s",LETTERS[i])
    keywords[['$TIMESTEP']]<-"0.0001"
    ##
    error.type<-match.arg(error.type)
    if(error.type=='extra keyword'){
      set.seed(i)
      keywords[[paste0(sample(LETTERS,4),collapse = "")]]<-paste0(sample(letters,4),collapse = "")
    }
    if(error.type=='$PN conflict'){
      detector.names.mod<-c("Detector a-A","Detector A-A","detector a-A")
      if(i %in% 1:3){
        colnames(mat)[grep("detector A-A",colnames(mat))]<-detector.names.mod[i]
      }
    }
    ##
    parms.adf<-data.frame(
      name = colnames(mat),
      desc = NA,
      minRange = apply(mat,2,min),
      maxRange = apply(mat,2,max)
    );parms.adf[['range']]<-parms.adf$maxRange-parms.adf$minRange
    parms.adf<-parms.adf[,c('name',"desc","range","minRange","maxRange")]
    rownames(parms.adf)<-paste0("$P",seq(nrow(parms.adf)))
    parms.adf$desc[!grepl("Time|[FS]SC",parms.adf$name)]<-paste0("stain",seq(parameters))
    ##
    if(error.type=='$PS conflict'){
      if(i==1){
        parms.adf[grepl("stain1$",parms.adf$desc),'desc']<-"Stain 1"
      }
    }
    ##
    volts.names<-paste0(rownames(parms.adf)[!grepl("Time",parms.adf$name)],"V")
    volts.vec<-stats::setNames(nm=volts.names,{set.seed(1337);sample(0:500,length(volts.names))})
    TYPE.names<-paste0(rownames(parms.adf),"TYPE")
    TYPE.vec<-stats::setNames(nm=TYPE.names,c("Time","Forward_Scatter","Side_Scatter",rep("Unmixed_Fluorescence",26)))
    keywords<-c(keywords,as.list(volts.vec),as.list(TYPE.vec))
    ##
    fcs<-methods::new(
      "flowFrame",
      exprs = mat,
      parameters = Biobase::AnnotatedDataFrame(parms.adf),
      description = keywords
    )
    ##
    file.out<-file.path(tempdir(),keywords[['$FIL']])
    flowCore::write.FCS(fcs,file.out)
  }
}

raw.spectral.synthetic.fcs.rand<-function(fname="CD1337 FakeFluor (Beads).fcs",n.minor=5,n.blips=5){
  mat.scatter<-matrix(
    data = c(
      stats::runif(n=4000*2,min=9.8E5,max = 1.2E6),
      stats::runif(200*2,min=0,max=4E6),
      stats::runif(800*2,min=(9.8E5*2),max=(1.2E6*2))
    ),
    nrow = 5000,
    byrow = TRUE,
    ncol = 2,
    dimnames = list(NULL,c('FSC-A','SSC-A'))
  )

  detectors<-c("UV"=16,"V"=16,"B"=14,"YG"=10,"R"=8)
  detectors<-unlist(sapply(names(detectors),function(i){paste0(i,seq(detectors[[i]]),"-A")}))

  mat.detectors<-matrix(
    data = stats::runif(n=5000*64,min=-500,max = 1500),
    nrow = 5000,
    ncol = 64,
    dimnames = list(NULL,detectors)
  )

  j.major<-sample(detectors,1)
  i.major<-sample(nrow(mat.detectors),2000,replace = FALSE)
  mat.detectors[i.major,j.major]<-stats::runif(n=2000,min=3E5,max = 4E5)

  j.minor<-sample(detectors[!detectors %in% j.major],n.minor)
  for(j in j.minor){
    mat.detectors[i.major,j]<-stats::runif(
      n=2000,
      min=sample(seq(1E4,1E5,by=1E4),1),
      max=sample(seq(1.25E5,2.75E5,by=1E4),1)
    )
  }
  j.blips<-sample(detectors[!detectors %in% c(j.major,j.minor)],n.blips)
  for(j in j.blips){
    mat.detectors[i.major,j]<-stats::runif(
      n=2000,
      min=sample(seq(1E3,1E4,by=1E3),1),
      max=sample(seq(1.25E4,2.75E4,by=1E3),1)
    )
  }

  mat<-cbind(mat.scatter,mat.detectors)

  parms.adf<-data.frame(
    name = colnames(mat),
    desc = NA,
    minRange = apply(mat,2,min),
    maxRange = apply(mat,2,max)
  );parms.adf[['range']]<-parms.adf$maxRange-parms.adf$minRange
  parms.adf<-parms.adf[,c('name',"desc","range","minRange","maxRange")]
  rownames(parms.adf)<-paste0("$P",seq(nrow(parms.adf)))

  keywords<-list('$FIL'=fname)
  keywords[['$CYT']]<-"InSilicometer"
  keywords[['$CYTSN']]<-"1337"
  keywords[['$PROJ']]<-"spectr_package"
  keywords[['$TUBENAME']]<-sub(".fcs","",fname)

  type<-paste0(rownames(parms.adf)[!grepl("[FS]SC", parms.adf$name)],"TYPE")
  keywords<-c(keywords,as.list(stats::setNames(rep("Raw_Fluorescence",length(type)),nm=type)))
  keywords[[paste0(rownames(parms.adf)[parms.adf$name=='FSC-A'],"TYPE")]]<-"Forward_Scatter"
  keywords[[paste0(rownames(parms.adf)[parms.adf$name=='SSC-A'],"TYPE")]]<-"Side_Scatter"

  ##
  fcs<-methods::new(
    "flowFrame",
    exprs = mat,
    parameters = Biobase::AnnotatedDataFrame(parms.adf),
    # parameters = methods::as(parms.adf, "AnnotatedDataFrame"),
    description = keywords
  )
  ##
  file.out<-file.path(tempdir(),keywords[['$FIL']])
  flowCore::write.FCS(fcs,file.out)
}
