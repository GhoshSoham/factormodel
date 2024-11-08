###
rfishbing.rj<-function(evdA,c,x){

      E<-evdA$vec
      l<-evdA$val

      y<-abs(t(E)%*%x)
      s<-sign(t(E)%*%x)
      d<-t(E)%*%c

      for(i in sample(1:length(x)) ) {
        q<-y^2/(1-y[i]^2)
        th<-rtheta( .5*(dim(E)[1]-3),(l[i]-t(q[-i])%*%l[-i]),
                        sum(s[-i]*(q[-i]^.5)*d[-i]),d[i] )
        y[i]<-sqrt(th)
        y[-i]<-sqrt((1-th)*q[-i])
        s[i]<-(-1)^rbinom(1,1,1/(1+exp(2*sqrt(th)*d[i])) )

                                     }
      x<-E%*%(y*s)
      x/sqrt(sum(x^2))
                              }

#####

rfishbing<-rfishbing.rj

#####
rtheta<-function(k,a,b,c) {

  d1<-1/2 ; d2<- 1+min(k, max( k-a,-1/2) )  
  lmx<-0
  if(a> k - d2 +1 ) { lmx<-  (k-d2+1)*(log(k-d2+1)-log(a))+ a-(k-d2+1) }
  
  bb<-max(0,b)
  thmx<-  c^2/(c^2+bb^2)  ; if( abs(c)==0 & bb==0) { thmx<-1 }
  
  lfg.mx  <- lmx +  
             (b*sqrt(1-thmx)+abs(c)*sqrt(thmx)) +
             log(1+exp(-2*abs(c)))

  cmet<-FALSE
  while(!cmet) {

    th<-rbeta(1,d1,d2)
    lfg.th <- (k-d2+1)*log(1-th)+a*th + 
              (b*sqrt(1-th)+abs(c)*sqrt(th)) + 
              log(1+exp(-2*sqrt(th)*abs(c))) 

    cmet<-( log(runif(1)) < lfg.th - lfg.mx )

                 }
  th                          }
#####


#####
ldbing.o2<-function(Z,A,B,C){ sum( diag( t(C)%*%Z + B%*%t(Z)%*%A%*%Z ) ) }

theta<-seq(0,2*pi,length=100)
O2<-array(dim=c(2,2,2*length(theta)))
O2[1,1,]<-c(cos(theta),cos(theta))
O2[2,1,]<-c(sin(theta),sin(theta))
O2[1,2,]<-c(sin(theta),-sin(theta))
O2[2,2,]<-c(-cos(theta),cos(theta))



rfishbing.O2<-function(A,B,C) {
  lpz<- apply(O2,3,ldbing.o2,A=A,B=B,C=C)
  O2[,,sample(1:length(lpz),1,prob=exp(lpz-max(lpz))) ]
                           }
#####
# X is a kxd matrix
d = 2
k = 20
A = matrix(0, k, k)
C = matrix(rnorm(4), k, d)
B = matrix(0, d, d)

sam = rfishbing.O2(A, C, B)



### kume and walker method
rs_gibbs<-function(s,l) {

  v<-runif(1, 0, exp(-sum(l*s)) )
  w<-runif(1, 0, (1-sum(s))^(-.5) )


  for(i in sample(1:length(s))) {

    c<- max( 0, 1-w^(-2)-sum(s[-i]) )
    d<- min( -(log(v)+sum(l[-i]*s[-i]))/l[i] , 1-sum(s[-i]))

    u<-runif(1, c^.5, d^.5)
    s[i]<-u^2
                                }
 
    s                           }


rfishbing.kw<-function(evdA,x) {

      E<- evdA$vec
      l<- -evdA$val

      l<-l-min(l)
      y<- t(E)%*%x
      s<- y^2
     
      tmp<- rs_gibbs(s[l>0],l[l>0])
      s[l==0]<- 1-sum(tmp)
      s[l>0]<-tmp
   
      y<-sqrt(s)*((-1)^rbinom(length(s),1,.5))
      x<-E%*%y
      x/sqrt(sum(x^2))           }
 




