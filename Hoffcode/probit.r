rzq_fc<-function(Z,pp_zq=1/100){
  sd_zq<-1/sqrt(pp_zq)
  zq<-c(-Inf,rep(NA,max(Ranks,na.rm=T)-1),Inf)
   for(ry in 1:(max(Ranks,na.rm=T)-1)){
    ub<-suppressWarnings(min(Z[ Ranks==ry+1 ],na.rm=T ) )
    lb<-suppressWarnings(max(Z[ Ranks==ry ],na.rm=T ) )
    zq[ry+1]<-  qnorm( runif(1,pnorm(lb,0,sd_zq),pnorm(ub,0,sd_zq)),0,sd_zq  )


    if( zq[ry+1]==Inf) { zq[ry+1]<-lb }
    if( zq[ry+1]==-Inf) { zq[ry+1]<-ub }
  
  
  
                                     }
 zq                          }
#####

#####
rZ_fc<-function(EZ,zq){

  for(ry in uRanks){
    ir<- ( Ranks==ry & !is.na(Ranks) )
    lb<- zq[ry]
    ub<- zq[ry+1]

    z<-qnorm(runif(sum(ir),pnorm(lb,EZ[ir],1),
                           pnorm(ub,EZ[ir],1)),EZ[ir],1)
    z[z== Inf]<-lb
    z[z==-Inf]<-ub
    z<-pmin( ub,pmax(lb,z))
    Z[ir]<-z
                           }
    ir<-is.na(Ranks)
    Z[ir]<-rnorm(sum(ir),EZ[ir],1)

    Z                            }
#####



#####
rZ.udm_fc<-function(EZ,zq,sig=1){
  
  Z<-EZ*0
  ut<-upper.tri(EZ)

  for(ry in uRanks){
    ir<- ( Ranks==ry & !is.na(Ranks) & ut )
    lb<- zq[ry]
    ub<- zq[ry+1]

    z<-qnorm(runif(sum(ir),pnorm(lb,EZ[ir],sig),
                           pnorm(ub,EZ[ir],sig)),EZ[ir],sig)
    z[z== Inf]<-lb
    z[z==-Inf]<-ub
    z<-pmin( ub,pmax(lb,z))
    Z[ir]<-z
                           }
    ir<-is.na(Ranks &ut )
    Z[ir]<-rnorm(sum(ir),EZ[ir],1)
    Z<-Z+t(Z)
#    diag(Z)<-rnorm(dim(Z)[1],0,sig*sqrt(2))
     diag(Z)<-rnorm(dim(Z)[1],diag(EZ),sig*sqrt(2))


    Z                            }
#####




addlines<-function(u1,u2,Y,col="green",lwd=1,lty=1) {
n<-dim(Y)[1]
for(i in 1:(n-1)){
for(j in (i+1):n){
if(!is.na(Y[i,j])) {
if(Y[i,j]!=0) {   segments(u1[i],u2[i],u1[j],u2[j],col=col,lwd=lwd) }
                     }
               }}
                              }

