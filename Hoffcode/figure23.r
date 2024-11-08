load("image.R3.rs0.m270")
ULU.pm0<-ULU.ps/pbs
LS0<-LS


load("image.R3.rs1.m270")
ULU.pm1<-ULU.ps/pbs
LS1<-LS



#### eigenvalues
pdf("eigenmix.pdf",height=4,width=8,family="Times")
par(mfrow=c(1,1),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
R<-dim(LS0)[2]
scans<-(1:dim(LS)[1])*100
plot( c(1,max(scans)),range(c(LS0,LS1)),type="n",xlab="iteration",ylab=
      expression(lambda))
for(r in 1:R){ lines((1:dim(LS1)[1])*100, LS1[,r],lwd=1,col="black") }
for(r in 1:R){ lines((1:dim(LS0)[1])*100, LS0[,r],lwd=1,col="gray") }

L.pm<- .5*(  sort(apply(LS0[-(1:100),],2,mean)) +  
             sort(apply(LS1[-(1:100),],2,mean))  )
abline(h=L.pm) ; abline(h=0,col="gray",lty=2)
dev.off()
####


####
YMLV<-dget("Butland_MLV")
rownames(Y)<- (dimnames(YMLV)[[1]])[match(colnames(Y),dimnames(YMLV)[[2]])]

pdf("upm.pdf",height=4,width=8,family="Times")
par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.75,.75,0))
ULU.pm<- (ULU.pm0+ULU.pm1)/2
tmp<-eigen(ULU.pm)
U<-tmp$vec[ ,c(1,2,m) ]
L<-diag(tmp$val[c(1,2,m)])

plot(U[,1:2]*1.1,type="n",xlab=expression(italic(U[1])),ylab=expression(italic(U[2])))
abline(v=0, col="black",lty=2) ; abline(h=0,col="black",lty=2)
addlines(U[,1],U[,2],Y,col="gray")
m12<-apply((U[,1:2])^2,1,sum)
b12<- (1:m)[ m12>=quantile(m12,.9) ]
s12<- (1:m)[ m12<quantile(m12,.9) ]
text(U[b12,1],U[b12,2], rownames(Y)[b12],cex=.8)
points(U[s12,1:2],pch=16,cex=.7)
bp3<-(1:m)[ U[,3]>=quantile(U[,3],.95) ]
bn3<-(1:m)[ U[,3]<=quantile(U[,3],.05) ]
b3<-sort(c(bp3,bn3))
points(U[b3,1:2],pch=21,bg="white")
text(U[b12,1],U[b12,2], rownames(Y)[b12],cex=.8)
####


####
s3<-(1:m)[-b3]
x<- apply(U[,1:2]%*%L[1:2,1:2],1,sum)
y<- U[,3]*L[3,3]
plot(x,y,type="n",xlab=expression(italic(lambda[1]*U[1]+lambda[2]*U[2]))  , 
      ylab=expression(italic(U[3])) )
addlines(x,y,Y,col="gray")
abline(h=0)
points(x[s3],y[s3],pch=16,cex=.7)
text(x[b3],y[b3], rownames(Y)[b3],cex=.8)
####

dev.off()



