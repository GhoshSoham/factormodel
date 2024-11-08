source("rfishbing.r")

######
rbing.kw<-function(x,lambda) 
{ 
  l<- -lambda
  l<-l-min(l)
  s<-x^2
  tmp<- rs_gibbs(s[l>0],l[l>0])
  s[l==0]<- 1-sum(tmp)
  s[l>0]<-tmp
  x<-sqrt(s)
  (x/sqrt(sum(x^2)))*(-1)^rbinom(length(lambda),1,.5) 
}
#######

#######
rbing.ph<-function(x,lambda) 
{ 
  for(i in sample(1:length(x)) ) 
  {
    q<-x^2/(1-x[i]^2)
    th<-rtheta( .5*(length(x)-3),(lambda[i]-t(q[-i])%*%lambda[-i]),0,0)
    x[i]<-sqrt(th) 
    x[-i]<-sqrt( (1-th)*q[-i] )
  }  
  (x/sqrt(sum(x^2)))*(-1)^rbinom(length(lambda),1,.5) 
}
######


######
rbing.sl<-function(x,lambda)
{
  for(i in sample(1:length(x)) )
  {
    q<-x^2/(1-x[i]^2)
    k<- .5*(length(x)-3) ; a<-lambda[i]-t(q[-i])%*%lambda[-i]
    u<-runif(1,0,exp(a*x[i]^2))

    if(a>0) {  lb<-pbeta( max(0, log(u)/a)  ,1/2,k+1) ; ub<-1 }
    if(a<0) {  ub<-pbeta( min(1, log(u)/a),  1/2,k+1) ; lb<-0 }

    th<-qbeta( runif(1,lb,ub) , 1/2,k+1  )
    x[i]<-sqrt(th)
    x[-i]<-sqrt( (1-th)*q[-i] ) 
  }
  (x/sqrt(sum(x^2)))*(-1)^rbinom(length(lambda),1,.5)
}
#######



########
rX.ph<-function(lambda,S) 
{
  p<-length(lambda) 
  X.ph<-matrix( rep(1/sqrt(p),p),S,p)
  set.seed(1)
  for(s in 2:S)
  {
    X.ph[s,]<- rbing.ph(X.ph[s-1,],lambda )
   }
  X.ph
} 
#####

#####
rX.sl<-function(lambda,S)
{
  p<-length(lambda)
  X.sl<-matrix( rep(1/sqrt(p),p),S,p)
  set.seed(1)
  for(s in 2:S)
  {
    X.sl[s,]<- rbing.sl(X.sl[s-1,],lambda )
   }
  X.sl
}
#####



#####
rX.kw<-function(lambda,S) 
{
  p<-length(lambda) 
  X.kw<-matrix( rep(1/sqrt(p),p),S,p)
  set.seed(1)
  for(s in 2:S)
  {
    X.kw[s,]<- rbing.kw(X.kw[s-1,],lambda )
   }
  X.kw
}
#####


library(coda)
RES<-NULL
for(k in c(2,4,8)) { 
for(g in c(0,5,10)) {
p<-k*10
lambda<-(p:1)/10 + g*(rep(c(1,0),c(k,p-k)) )

S<-5000

t.kw<-system.time( X.kw<-rX.kw(lambda,S) )[3]
t.ph<-system.time( X.ph<-rX.ph(lambda,S) )[3]

s.kw<-mean(apply(X.kw[,1:5]^2,2,effectiveSize))
s.ph<-mean(apply(X.ph[,1:5]^2,2,effectiveSize))


tmp<-c(p,g,s.kw,s.ph,s.kw/t.kw,s.ph/t.ph,s.ph*t.kw/(s.kw*t.ph))
RES<-rbind(RES,tmp)
cat(round(tmp,2),"\n") 
                        } }

tmp<-  cbind(  RES[,4]/RES[,3] , RES[,7] )
tmp<-cbind( tmp[1:3,], tmp[4:6,],tmp[7:9,] )
print( round(tmp,2))


