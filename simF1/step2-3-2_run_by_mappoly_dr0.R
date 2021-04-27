# install.packages("devtools")
# devtools::install_github("jendelman/diaQTL")
# devtools::install_github("mmollina/mappoly", dependencies=TRUE)
library(mappoly)
packageVersion("mappoly")
packageVersion("diaQTL")

library(rstudioapi)    
workdir0=dirname(rstudioapi::getActiveDocumentContext()$path)
workdir=file.path(workdir0,"res_mappoly")
setwd(workdir)
getwd()
source("../mappoly2polyorigin_V2.R")

# run mappoly
dr = 0
sumfile = paste0("ressum_F1_mappoly_DR",dr,".txt")
res <- c("mappoly_version",paste0(packageVersion("mappoly")))
write(res, file = sumfile,ncolumns = length(res),append = TRUE, sep = ",")
telap=c("nfilt","tread","tmakeseq","tpair","thmm","thmmerr","tprob","tsave")
res<-c(c("popsize","DR","rep"),telap)
write(res, file = sumfile,ncolumns = length(res),append = TRUE, sep = ",")
nrep <-50
for (dr in c(dr)){
  # see https://github.com/mmollina/Test_mappoly/ for mrk.tail and n.ph
  mrk.tail <- c(60, 60, 60, 100, 100, 150, 200)
  n.ph <- c(20, 20, 40, 60, 60, 100, 300)
  popsizels=c(200,100,50,30,20,15,10)
  for (i in 1:6){
    popsize = popsizels[i]
    for (rep in 1:nrep){
      dataid=paste0("F1_noff",popsize,"_DR",dr,"_rep",rep,"_")
      print(dataid)
      genofile = paste0(dataid,"polyorigin_geno_snparray_mappoly.csv")
      outfile = paste0(dataid,"genoprob_mappoly.csv")
      tread=system.time(dat.dose.csv <- read_geno_csv(file.in=genofile, ploidy=4))
      tmakeseq=system.time({
        pval.bonf <- 0.05/dat.dose.csv$n.mrk
        dat.chi.filt <- filter_segregation(dat.dose.csv, chisq.pval.thres =  pval.bonf, 
                                           inter = FALSE)
        nfilt=length(dat.chi.filt$exclude)
        dat.seq <- make_seq_mappoly(dat.chi.filt) 
      })
      tpair=system.time(all.rf.pairwise <-  est_pairwise_rf(input.seq = dat.seq))
      eps0 <- 0.01
      thmm=system.time(
        map <- est_rf_hmm_sequential(input.seq = dat.seq,
                                     start.set = 4,
                                     thres.twopt = 5,
                                     thres.hmm = 5,
                                     extend.tail = mrk.tail[i],
                                     twopt = all.rf.pairwise,
                                     sub.map.size.diff.limit = 10,
                                     phase.number.limit = n.ph[i]))
      thmmerr=system.time(map.error <- est_full_hmm_with_global_error(input.map = map, error = eps0))
      tprob=system.time(genoprob <- calc_genoprob_error(input.map = map.error, error = eps0))
      tsave=system.time(mappoly2polyorigin(map.error,genoprob,outfile=outfile))
      telap <- as.numeric(c(nfilt,tread[3],tmakeseq[3],tpair[3],thmm[3],thmmerr[3],tprob[3],tsave[3]))
      res<-c(c(popsize,dr,rep),telap)
      print(res)
      write(res, file = sumfile,ncolumns = length(res),append = TRUE, sep = ",")
    }
  }
}  
