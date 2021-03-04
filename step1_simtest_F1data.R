# Installing and using the PedigreeSimR to simulate the cross and GBS data
# devtools::install_github("rramadeu/PedigreeSimR")
library(PedigreeSimR)
packageVersion("PedigreeSimR")

library(rstudioapi)    
workdir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
getwd()
df = read.table("TableS1_RussetFounderHaplo.csv",sep=",",header = TRUE,skip=1)
setwd("F1data")
russet.map = df$cM
russet.hap = matrix(unlist(df[,-c(1:2)]), nrow=length(russet.map))-1
ploidy=4
np=2
eps0=0.01
for (popsize in c(10,15,20,30,50,100,200))
  for (dr in c(0.0, 0.5))
    for (rep in 1:50){
      outid=paste0("F1_noff", popsize,"_DR",dr,"_rep",rep,"_")
      pedigree = diallel_pedigree(parents=np,popsize=popsize)
      pedigreesimR(russet.map,russet.hap[,1:(np*ploidy)],sampleHap = FALSE,filename=outid,
                   pedigree,ploidy=ploidy,
                   prefPairing = dr,
                   quadrivalents = dr,
                   workingfolder = getwd(),
                   epsilon = c(eps0,eps0),
                   missingFreq=c(0,0.1),
                   trackErrorSim=FALSE
      )
      sapply(list.files(pattern = paste0(outid,"pedsim*")), unlink)
      unlink("*truevalue.csv")
      unlink("*geno.csv")
    }


