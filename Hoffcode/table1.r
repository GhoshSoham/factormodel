source("rlangevin.r")

RES<-NULL
##
for(m in c(10,20,30)) {
for(R in  c(2,4,6)) {
for(d in  m*c(.5,1,2)) {
  set.seed(1)
  svdM<-svd(matrix(rnorm(m*R),m,R))
  M<-svdM$u[,1:R]%*%diag(rep(d,R))%*%t( svdM$v[,1:R] ) 
  nrej<-NULL    
  for(s  in 1:100) { 
  nrej<-c(nrej, rlangevin(M)$rej ) 
#  cat(m,R,d,s,mean(nrej),"\n")
                    }
  res<-c(m,R,d,mean(nrej))
  RES<-rbind(RES,res) 
  cat(res,"\n")          }}}

###

plot(c(20,100),range(log(RES[,4])),type="n")
for( R in c(2,4,6)) {
tmp<-RES[RES[,2]==R ,]
lines(tmp[,3],log(tmp[,4]) )
                     }

###




