library(rstudioapi)    
workdir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
getwd()

# prepare input genofile for mappoly
popsizels=c(10,15,20,30,50,100,200)
for (dr in c(0,0.5))
  for (popsize in popsizels)
    for (rep in 1:50){
      dataid=paste0("F1_noff", popsize,"_DR",dr,"_rep",rep,"_")
      genofile = paste0("F1data/",dataid,"polyorigin_geno_snparray.csv")
      tmp = read.csv(genofile,strip.white=TRUE,sep=',')
      tmp = tmp[,c(1,4,5,2,3,6:ncol(tmp))]
      genofile2 = paste0("res_mappoly/",dataid,"polyorigin_geno_snparray_mappoly.csv")
      write.csv(tmp,file=genofile2,row.names=FALSE)
    }