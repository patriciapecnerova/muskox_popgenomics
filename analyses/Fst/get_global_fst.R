##define Fst function
getFst<-function(est){
  N1<-nrow(est)-1
  N2<-ncol(est)-1
  cat("N1: ",N1 ," N2: ",N2,"\n")
  est0<-est
  est0[1,1]<-0
  est0[N1+1,N2+1]<-0
  est0<-est0/sum(est0)
  
  aMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
  aMat.ss<<-matrix(NA,nrow=N1+1,ncol=N2+1)
  baMat<<-matrix(NA,nrow=N1+1,ncol=N2+1)
  for(a1 in 0:(N1))
    for(a2 in 0:(N2)){
      p1 <- a1/N1
      p2 <- a2/N2
      q1 <- 1 - p1
      q2 <- 1 - p2
      N <- (p1-p2)^2
      D <- p1*(1-p2)+p2*(1-p1)
      aMat[a1+1,a2+1]<<-N
      baMat[a1+1,a2+1]<<-D
      #sample size correction
      N.ss <- (p1-p2)^2-((p1*(1-p1))/(N1-1))-((p2*(1-p2))/(N2-1))
      aMat.ss[a1+1,a2+1]<<-N.ss
    }
  ## sample size corrected moment estimator
  ss <- sum(est0*aMat.ss,na.rm=T)/sum(est0*baMat,na.rm=T)
  c(fstSS=ss)
}

## insert population names instead of "pop1" and "pop2" for all pairs
full<-scan("muskox_bamlist_pop1.pop2.sfs")
#N of individuals in the group*2+1
fullM<-matrix(full, nrow=55, ncol=31, byrow=T)
print(getFst(fullM))
