# library(diaQTL)

mappoly2polyorigin<-function(input.map,genoprob,outfile="genoprob_mappoly.csv")
{
  data(F1codes,package="diaQTL") ##Line Changed
  ## Exports parental phasing
  p = input.map$maps[[1]]$seq.ph[[1]]
  q = input.map$maps[[1]]$seq.ph[[2]]
  inds=dim(genoprob$probs)[3]
  markers=dim(genoprob$probs)[2]
  PQmat = matrix(0,8,length(p))
  colnames(PQmat) = dimnames(genoprob$probs)[[2]]
  pq.out = matrix(NA,length(p),2)
  for(i in 1:length(p)){
    if(p[[i]][1]>0){
      PQmat[1:4,i][p[[i]]] = 1
      pq.out[i,1] = paste0(PQmat[1:4,i]+1,collapse="|")
    }else{
      pq.out[i,1] = "1|1|1|1"
    }
    if(q[[i]][1]>0){
      PQmat[5:8,i][q[[i]]] = 1
      pq.out[i,2] = paste0(PQmat[5:8,i]+1,collapse="|")
    }else{
      pq.out[i,2] = "1|1|1|1"
    }
  }
  pq.out <- cbind(dimnames(genoprob$probs)[[2]],"A",round(genoprob$map,4),pq.out)
  colnames(pq.out) <- c("marker", "chromosome" , "position" , "A" , "B")
  # write.csv(pq.out,file= outfile,quote=FALSE,row.names=FALSE)
    ## Exports offspring phasing
  F1states = F1codes$State
  F1states = stringr::str_split(F1states,"-",simplify=TRUE)
  F1states = apply(F1states,2,as.numeric)
  F1states = matrix(letters[F1states],ncol=4)
  F1states = paste0(F1states[,1],F1states[,2],":",F1states[,3],F1states[,4])
  mappolycodex = match(rownames(genoprob$probs),F1states)
  options(scipen=999)
  genoprob$probs = round(genoprob$probs,4)
  OS.out = matrix(NA,markers,inds)
  colnames(OS.out) =  dimnames(genoprob$probs)[[3]]
  for(i in 1:inds){
    for(j in 1:markers){
      keep = which(genoprob$probs[,j,i]>0)
      OS.out[j,i] =  paste0(paste0(mappolycodex[keep],collapse="|"),
                            "=>",
                            paste0(genoprob$probs[keep,j,i],collapse="|"))
    }
  }
  OS.out = cbind(pq.out,OS.out)
  write.csv(OS.out,file= outfile,quote=FALSE,row.names=FALSE)
}
