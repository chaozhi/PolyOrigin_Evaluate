# install.packages("devtools")
# devtools::install_github("jendelman/diaQTL")
# devtools::install_github("mmollina/mappoly", dependencies=TRUE)
library(mappoly)
packageVersion("mappoly")
packageVersion("diaQTL")

library(rstudioapi)    
workdir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
getwd()
source("mappoly2polyorigin_V2.R")

# run mappoly
for (popindex in 1:3){
  pop=paste0("pop",popindex)
  sumfile = paste0("ressum_mappoly_",pop,".txt")
  res <- c("mappoly_version",paste0(packageVersion("mappoly")))
  write(res, file = sumfile,ncolumns = length(res),append = TRUE, sep = ",")
  telap=c("nfilt","tread","tmakeseq","tpair","thmm","thmmerr","tprob","tsave")
  res<-c(c("dataid","chr","pop"),telap)
  write(res, file = sumfile,ncolumns = length(res),append = TRUE, sep = ",")
  for (chr in 1:12){
    dataid=paste0("potato_",pop, "_chr",chr)
    print(dataid)
    genofile = paste0("..//step1_data//",dataid,"_dose_mappoly.csv")
    # file.exists(genofile)
    outfile = paste0(dataid,"_genoprob_mappoly.csv")
    tread=system.time(dat.dose.csv <- read_geno_csv(file.in=genofile, ploidy=4))
    #here assumes that the input order is the final order (skips missing filtering and grouping)
    tmakeseq=system.time({
      pval.bonf <- 0.05/dat.dose.csv$n.mrk
      dat.chi.filt <- filter_segregation(dat.dose.csv, chisq.pval.thres =  pval.bonf, 
                                         inter = FALSE)
      nfilt=length(dat.chi.filt$exclude)
      dat.seq0 <- make_seq_mappoly(dat.chi.filt) 
    })
    tpair=system.time({
      tpt <-  est_pairwise_rf(input.seq = dat.seq0)
      dat.seq  <- rf_snp_filter(tpt)
    })
    thmm=system.time(map <- est_rf_hmm_sequential(input.seq = dat.seq,
                                                       start.set = 3,
                                                       thres.twopt = 10,
                                                       thres.hmm = 10,
                                                       extend.tail = 50,
                                                       twopt = tpt,
                                                       verbose = TRUE,
                                                       tol = 10e-2,
                                                       tol.final = 10e-3,
                                                       phase.number.limit = 40,
                                                       sub.map.size.diff.limit = 5,
                                                       info.tail = TRUE,
                                                       reestimate.single.ph.configuration = TRUE))
    map<-filter_map_at_hmm_thres(map, thres.hmm = 0.001)
    eps0=0.02
    thmmerr=system.time(map.error <- est_full_hmm_with_global_error(input.map = map, error = eps0))
    tprob=system.time(genoprob <- calc_genoprob_error(input.map = map.error, error = eps0))
    tsave=system.time(mappoly2polyorigin(map.error,genoprob,outfile=outfile))
    telap <- as.numeric(c(nfilt,tread[3],tmakeseq[3],tpair[3],thmm[3],thmmerr[3],tprob[3],tsave[3]))
    res<-c(c(dataid,chr,pop),telap)
    print(res)
    write(res, file = sumfile,ncolumns = length(res),append = TRUE, sep = ",")
  }
}