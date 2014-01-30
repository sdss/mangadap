#include <stdio.h>

void bvls(int argc, void *argv[]) {

extern void bvls_(); /* Fortran routine */
int *n, *m, *nsetp, *index, *ierr;     /* These are scalars */
double **a, **b, **bnd, **w, **x;   /* These are arrays */
double *rnorm;


/* This block makes sure that the variables
going into the F90 routine have the right type */
a =     (double **) argv[0];
m =     (int *)     argv[1];
n =     (int *)     argv[2];
b =     (double **) argv[3];
bnd =   (double **) argv[4];
x =     (double **) argv[5];
rnorm = (double *)  argv[6];
nsetp = (int *)     argv[7];
w =     (double **) argv[8];
index = (int *)     argv[9];
ierr =  (int *)     argv[10];

bvls_(a,m,n,b,bnd,x,rnorm,nsetp,w,index,ierr); /* Run the F90 BVLS */

}



