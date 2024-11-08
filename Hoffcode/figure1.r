Y<-dget("Y_Pro")

library(network) 

addlines<-function(u1,u2,Y,col="green",lwd=1,lty=1) {
n<-dim(Y)[1]
for(i in 1:(n-1)){
for(j in (i+1):n){
if(!is.na(Y[i,j])) {
if(Y[i,j]!=0) {   segments(u1[i],u2[i],u1[j],u2[j],col=col,lwd=lwd) }
                     }
               }}
                              }



xy<-plot(as.network(Y))
dg<-apply(Y,2,sum,na.rm=T)

pdf("prodesc.pdf",height=4,width=8,family="Times") 

par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))

plot(xy,type="n",xlab="",ylab="",xaxt="n",yaxt="n")
addlines(xy[,1],xy[,2],Y,col="gray") 
points(xy,pch=16,col="black",cex=.7)

plot( 1:max(dg), ( table( c(dg, 1:max(dg)))-1 ) ,type="h",lwd=3,
      xlab="degree",ylab="# nodes with a given degree" )


dev.off()

