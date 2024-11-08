#include<math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void rW(double *kap, int *m, double *W)
{

        GetRNGstate();


        double b=( -2.0*(*kap) + sqrt(  4*pow(*kap,2)+pow(*m-1.0,2) ) )/(*m-1.0) ;

        double x0=(1.0-b)/(1.0+b) ;
        double c= (*kap)*x0 +(*m-1.0)*log(1.0-pow(x0,2)) ;

        double Z, U;

        int done=0;
        while( done==0) {

                Z=rbeta( (*m-1.0)/2.0, (*m-1.0)/2.0 );
                *W= ( 1-(1+b)*Z)/(1.0-(1.0-b)*Z) ;
                U=runif(0,1);
                if(   (*kap)*(*W)+(*m-1)*log(1-x0*(*W))-c  > log(U) )  {done=1;}
                                 }
        PutRNGstate();

}

