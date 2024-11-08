source("rfishbing.r") ;  source("probit.r") ; library(MASS) ; library(abind)
Y<-dget("Y_Pro")


rs<-0
R<-3
m<-dim(Y)[1]
set.seed(1)
ranks<-rank(c(Y),ties.method="random",na.last="keep")
Z<-array( qnorm(ranks/(sum(!is.na(ranks))+1)),dim=dim(Y))
Z[lower.tri(Z)]<-0 ; Z<-Z+t(Z)  ;  diag(Z)<-diag(Z)/2
Z[is.na(Z)]<-rnorm(sum(is.na(Z)) )
uRanks<-1:length(unique(c(Y[!is.na(Y)])))
Ranks<-array(match(c(Y), sort(unique(c(Y)))),dim=dim(Y))
zq<-rzq_fc(Z)
tmp<- eigen( .5*(Z+t(Z)))
U<-tmp$vec[ ,order(-abs(tmp$val))[1:R] ]
L<-diag(tmp$val[order(-abs(tmp$val))[1:R] ])
s2<-1
t2<-m
seed<-1
if(rs==1) {
  set.seed(seed)
  L<-L*0
  U<-svd(matrix(rnorm(m*m),m,m))$u[,1:R] 
          }
#########
LS<-matrix(nrow=0,ncol=R)
pbs<-0
ULU.ps<-matrix(0,m,m)
BURN<-10000
NSCAN<-110000
for(s in 1:NSCAN) {

  zq<-rzq_fc(Z)
  Z<-rZ.udm_fc(U%*%L%*%t(U),zq)

  for(r in 1:R) { 
    L[r,r]<-rnorm(1,t(U[,r])%*%Z%*%U[,r]/(2*s2/t2+1),1/sqrt(1/t2 +1/(2*s2) ))
                }

  for(r in 1:R) {    
    N<-Null(U[,-r])
    M<- t(N)%*%(Z/2)%*%N*L[r,r]/s2
    U[,r]<-N%*%rfishbing(eigen(M),rep(0,m-R+1),t(N)%*%U[,r] )
                 }

  #### 
  if(s%%100==0) {
    cat(s,t2," ",round(diag(L),1),diag(t(U)%*%U),"\n")
    LS<-rbind(LS,diag(L))
    #plot(c(1,dim(LS)[1]),range(LS),type="n") ;apply(LS,2,lines) 
    #abline( h=apply(LS,2,mean)) ; abline(h=0,col="gray",lty=2)
    
    if(s>BURN) { ULU.ps<-ULU.ps+U%*%L%*%t(U) ; pbs<-pbs+1 }
    if(s%%1000==0) { save.image(paste("image.R",R,".rs",rs,".m",m,sep=""))  }
               }
   ####
    
                  }

#####




