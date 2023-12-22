#include <stdio.h>  //printf
#include <stdlib.h> //malloc
#include <math.h>   //sqrt


//#include <lapacke.h> //linear algebra for Newton-Raphson
#include <gsl/gsl_cblas.h> //cblas_ routines

#include "clsqr2/lsqr.h" //lsqr for linear algebra

//Macros
#define VECNORM(v) sqrt( v[0]*v[0] + v[1]*v[1] + v[2]*v[2] )
#define DOT3D(u,v) ( u[0]*v[0] + u[1]*v[1] + u[2]*v[2] )
#define HELPER_JACOBIAN(r,d_c,d_p,i,j) (-(i==j)/d_c + 3*r[i]*r[j]/d_p) //component of Jacobian of 1/r^2 force
#define IDX(i,j,m,n)  (i)*(n)+(j)

//Symplectic integration macros for readability
#define UPDATE_Q  cblas_daxpy( 6, c[stage], p, 1, q, 1 )
#define UPDATE_P  cblas_daxpy( 6, d[stage], f, 1, p, 1 )
#define UPDATE_DQ if( dx != NULL){ cblas_daxpy( 6*6, c[stage], C, 1, A, 1 ); cblas_daxpy( 6*6, c[stage], D, 1, B, 1 ); }
#define UPDATE_DP if( dx != NULL){  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 6, 6, d[stage], df, 6, A, 6, 1.0, C, 6); cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 6, 6, 6, d[stage], df, 6, B, 6, 1.0, D, 6); }




//Store parameters in a struct for easy access
typedef struct par{

  //Integration parameters
  double min_r; //Don't simulate anything below this radius
  double max_r; //Don't simulate anything above this radius
  int N; //number of timesteps

  //Newton-Raphson parameters
  int maxit;    // maximum Newton iterations
  double rcond; // 
  double damp;  // weight Newton steps by this number
  double acceptable_F;
}par;






//function prototypes
int newton_raphson( double *z, par *params );
int symplectic_integration( double *x, double *dx, double dt, par *p );
int force( double *f, double *df, double *x, par *p );
int objective( double *F, double *dF, double *z, par *params);
double hamiltonian( double *x );
int rescale( double *x );

void print_vector( double *x, int N );
void print_matrix( double *x, int m, int n );



int main( int argc, char *argv[] ){
  double z[13] = {  1,0,0,   //r1
                    0,0,0.1, //r2
                    0,1,0,   //p1
                    0,0,0,   //p2
                    7}; //T
 
  double dx[12*12];
  double dF[12*13];
  double F[12];


  //Set up parameters
  par params;
  params.min_r = 0.01; //Halt integration if distance becomes this small
  params.max_r = 10.0; //Halt integration if distance becomes this big
  
  params.N     = pow(2,10);
  params.maxit = pow(2,11);
  params.rcond = -1; //deprecated parameter
  params.damp  = 1e-2;
  params.acceptable_F = 1e-3;

  int flag;
  for( int i=0; i<1024*1024; i++){
    //Pick a random state
    for(int j=0; j<12; j++){
      z[j] =  2*( (float) rand() / RAND_MAX) - 1; // between [-1,1)
    }
    //Rescale state such that H=-1 before setting period
    if( rescale( z ) != 0 ){ continue; } //pick a new state if not a bound stae
    z[12] = 15*( (float) rand() / RAND_MAX); // pick a period between [0,15)

    flag = newton_raphson( z, &params );
 
    if(flag < 0) continue; //solution not found

    char filename[32];
    sprintf( filename, "states/%08d.txt", i);
    FILE *ptr = fopen( filename, "w" );

    for( int j=0; j<13; j++){
      fprintf( ptr, "%.16e\n", z[j] );
    }
    fclose(ptr);
  }
  return 0;
}








void print_vector( double *x, int n ){
  printf("v = \n" );
  for( int i=0; i<n; i++){
    printf("%.12f\n", x[i] );
  }
  printf("\n\n");
}


void print_matrix( double *x, int m, int n ){
  printf("M = \n" );
  for( int i=0; i<m; i++){
  for( int j=0; j<n; j++){
    printf("%.4f  ", x[IDX(i,j,m,n)] );
  }
  printf("\n");
  }
  printf("\n\n");
}


void aprod(int mode, int m, int n, double x[], double y[], void *UsrWrk ){
//     Function needed as an argument to lsqr
//     From documentation:
//     
//     The matrix A is intended to be large and sparse.  It is accessed
//     by means of subroutine calls of the form
//
//                aprod ( mode, m, n, x, y, UsrWrk )
//
//     which must perform the following functions:
//
//                If mode = 1, compute  y = y + A*x.
//                If mode = 2, compute  x = x + A(transpose)*y.
//
//     The vectors x and y are input parameters in both cases.
//     If  mode = 1,  y should be altered without changing x.
//     If  mode = 2,  x should be altered without changing y.
//     The parameter UsrWrk may be used for workspace as described
//     below.
//
//     I will use workspace to contain the matrix values.

  double *A = (double *) UsrWrk;
  
  if( mode == 1 ){
    // y = y + A*x
    for(int i=0; i<m; i++){
      for( int j=0; j<n; j++){
        y[i] += A[IDX(i,j,m,n)]*x[j];
      }
    }
  }
  
  
  if( mode == 2 ){
    // x = x + A(transpose)*y
    for(int i=0; i<m; i++){
      for( int j=0; j<n; j++){
        x[j] += A[IDX(i,j,m,n)]*y[i];
      }
    }
  }
}
 

int newton_raphson( double *z, par *params ){
  //Do Newton-Raphson to find a root of
  
  double z0[13]; //z will be overwritten on evaluation of F
  
  //Even though F is only 12D truly, the solution to the underdetermined system
  //will be stored in F. Therefore we need memory for 13 elements.
  double F[13];

  double *dF = (double *) malloc( sizeof(double)*12*13 );

  double F_norm;

  for( int i=0; i<(params->maxit); i++){
    cblas_dcopy(13, z, 1, z0, 1); //copy current value of z into z0;
  
    //Evaluate F and its Jacobian
    //Exit if there are integration errors
    if( objective( F, dF, z, params) != 0){ return -1; }
    F[12] = 0; //Manually set the final element of F to zero.

    //print the error at this stage to screen
    F_norm = cblas_dnrm2( 12, F, 1 );

    // Solve for the Newton step with the LSQR algorithm, which has a lightweight C implementation
    // See clsqr2/lsqr.c for documentation
    int m = 12;
    int n = 13;
    double damp = 0; //NOT the same damp as Newton step damp
    double *u = F; //RHS vector
    double v[13]; //workspace
    double w[13]; //workspace
    double x[13]; //returns computed solution! (Newton step in this case)

    double *se = NULL; //NULL to skip estimating standard errors


    //Documentation says set to zero if you don't particularly care
    double atol   = 1e-16; //relative error in matrix
    double btol   = 1e-16; //relative error in rh
    double conlim = 1e32;
    int itnlim    =   20; //iteration limit
    FILE *nout    = NULL; //stdout

    void * UsrWrk = (void *) dF; //pass in a poitner to matrix

    //Output memory for diagnostics
    int    istop_out;
    int    itn_out;
    double anorm_out;
    double acond_out;
    double rnorm_out;
    double arnorm_out;
    double xnorm_out;


    //Call lsqr
    lsqr( m, n, (*aprod), damp, UsrWrk,
      u,    // len = m
      v,    // len = n
      w,    // len = n
      x,    // len = n
      se,   // len = *
      atol,
      btol,
      conlim,
      itnlim,
      nout,
      // The remaining variables are output only.
      &istop_out,
      &itn_out,
      &anorm_out,
      &acond_out,
      &rnorm_out,
      &arnorm_out,
      &xnorm_out
     );

    //Print diagonostic info to screen
    printf("Newton step %d: |F| = %.6e, |r| = %.6e, |step| = %.6e \n", i, F_norm, rnorm_out, xnorm_out);
   
    if( F_norm != F_norm ){
      //we have nan. This is a trick that might not work for all compilers
      return -9;
    }
 
    // Update z. Keep in mind it was overwritten in objective(...);
    cblas_dcopy(13, z0, 1, z, 1);      // z <- z0
    cblas_daxpy(13, -(params->damp), x, 1, z, 1); // z <- z - damp*x, where x has been overwritten to be the Newton step

    //Rescale state after Newton step
    rescale( z );
  }

  //Cleanup memory
  free(dF);
  
  if( F_norm > (params->acceptable_F) ){ return -2; } //too high error
  if( fabs(z[12]) < 0.1 ){ return -3; } //too short of period

  //If we pass all checks, congratulations! You found a periodic orbit
  return 0;
}   


int symplectic_integration( double *x, double *dx, double dt, par *params ){
  /* PURPOSE:
   * This function integrates the three body problem with equal masses m=1 in 3D with a
   * symplectic integration scheme. Both the final position in state space and its Jacobian are computed.
   *
   * The special structure of symplectic integration is abused to compute the Jacobian very efficiently. 
   * BLAS routines are used for linear algebra. This probably makes no difference for adding vectors in 6D.
   *
   * The symplectic scheme was found by performing Newton's method on the thrid-order symplectic constraints. 
   * A third order method was converged with coefficients approximately given by
   *
   * c = { 0.2638, -0.2638, 1 };
   * d = { 0.8583, -0.1319, 0.2736 };
   *
   * The exact algebraic values are known and computed in the following code
   */

  double f[6];    //memory for force
  double df[6*6]; //memory for jacobian of force

  //pointers to position (q) and momentum (p)
  double *q = &x[0];
  double *p = &x[6];

  //Coefficient vectors for symplectic integration
  double c[3];
  double d[3];
  c[0] = (sqrt(7.0/3.0) - 1.0)/2;
  c[1] = -c[0];
  c[2] = 1.0;
  d[0] = (1+c[0])/(1-c[0])/2;
  d[1] = -c[0]/2;
  d[2] = 1 - d[0] - d[1];

  //Print out symplectic integration coefficients if you want
  //printf("%+0.5f\t%+0.5f\t%0.5f\n", c[0], c[1], c[2]);
  //printf("%+0.5f\t%+0.5f\t%0.5f\n", d[0], d[1], d[2]);

  //We might as well rescale c and d by the timestep
  for(int i=0; i<3; i++){
    c[i] *= dt;
    d[i] *= dt;
  }

  //Let the jacobian dx be [A, B; C, D] in MATLAB notation
  //We will eventually fill the Jacobian, but let's declare submatrices t.
  double *A;
  double *B;
  double *C;
  double *D;

  //initialize Jacobian
  if( dx != NULL){
    A = (double *) malloc( 6*6*sizeof(double) );
    B = (double *) malloc( 6*6*sizeof(double) );
    C = (double *) malloc( 6*6*sizeof(double) );
    D = (double *) malloc( 6*6*sizeof(double) );
 
    for(int i=0; i<6; i++){
      for( int j=0; j<6; j++){
        A[IDX(i,j,6,6)] = (i==j); // A = dq/dq
        B[IDX(i,j,6,6)] = 0;      // B = dq/dp
        C[IDX(i,j,6,6)] = 0;      // C = dp/dq
        D[IDX(i,j,6,6)] = (i==j); // D = dp/dp
      }
    }
  }


  //Time integration
  for( int i=0; i<(params->N); i++ ){
    for( int stage=0; stage<3; stage++){
      //I have hidden all of the details of this loop in UPDATE_ macros that call cblas routines.
      // q <- q + c[stage] * p
      UPDATE_Q;
      UPDATE_DQ;

      //compute the force and check for errors
      if( force(f, df, x, params) != 0 ){ return -1; } 
    
      // p <- p + d[stage] * f
      UPDATE_P;
      UPDATE_DP;
    }
  } //end timestepping
  

  //overwrite the Jacobian matrix
  if( dx != NULL){
    for(int i=0; i<6; i++){
      for(int j=0; j<6; j++){
          dx[IDX(i  ,j  ,12,12)] = A[IDX(i,j,6,6)];
          dx[IDX(i  ,j+6,12,12)] = B[IDX(i,j,6,6)];                                
          dx[IDX(i+6,j  ,12,12)] = C[IDX(i,j,6,6)];
          dx[IDX(i+6,j+6,12,12)] = D[IDX(i,j,6,6)];
      }
    }
  }

  //Free memory
  free(A);
  free(B);
  free(C);
  free(D);
}





int objective( double *F, double *dF, double *z, par *params){
  /* PURPOSE:
   * The objective function for periodic orbits: x(T) - x(0) = 0
   */
  
  //save initial state since z will be rewritten
  double x0[12]; //initial state
  cblas_dcopy( 12, z, 1, x0, 1 ); //copy initial state
  
  //shoot forward in time
  double  dt = z[12]/( params->N );
  double *dx = (double *) malloc( 12*12*sizeof(double) );
  if( symplectic_integration( z, dx, dt, params ) != 0){ return -1; }
  
  // F <- z[1:12] - x0 in matlab
  for(int i=0; i<12; i++){
    F[i] = z[i] - x0[i];
  }

  // dF = [ dx - I, v ]
  for(int i=0; i<12; i++){
    for(int j=0; j<12; j++){
      dF[IDX(i,j,12,13)] = dx[IDX(i,j,12,12)] - (i==j); //subtract identity
    }
  }
  
  double ff[6]; //memory for force
  force( ff, NULL, z, params ); //compute force. Pass NULL to skip Jacobian calculation
  
  for( int i=0; i<6; i++){
    dF[IDX(i,  12,12,13)] = z[i+6]; //time derivative of position is momentum
    dF[IDX(i+6,12,12,13)] = ff[i];  //time derivative of momentum is force
  }

  //free memory
  free(dx);

  return 0;
}




int force( double *f, double *df, double *x, par *p ){
    /*
    PURPOSE:
    Compute the force
    */

  //Separation vectors
  double r12[3] = {   x[0] -   x[3],   x[1] -   x[4],   x[2] -   x[5] }; //r1 - r2
  double r13[3] = { 2*x[0] +   x[3], 2*x[1] +   x[4], 2*x[2] +   x[5] }; //r1 - r3
  double r23[3] = {   x[0] + 2*x[3],   x[1] + 2*x[4],   x[2] + 2*x[5] }; //r2 - r3

  //Compute magnitudes
  double d12 = VECNORM(r12);
  double d13 = VECNORM(r13);
  double d23 = VECNORM(r23);

  //Find minimum d with two comparisons
  //If the minimum distance is too small, return an error
  double min;
  min = (d12 < d13) ? d12 : d13;
  min = (min < d23) ? min : d23;
  if( min < p->min_r ){
    //Failure: we have a very small separation in our trajectory
    return -1;
  }

  double max;
  max = (d12 > d13) ? d12 : d13;
  max = (max > d23) ? max : d23;
  if( max > p->max_r ){
    //Failure: we have a very large separation in our trajectory
    return -2;
  }



  //Cube distances
  double d12_c = d12*d12*d12;
  double d13_c = d13*d13*d13;
  double d23_c = d23*d23*d23;

  f[0] = -r12[0]/d12_c -r13[0]/d13_c;
  f[1] = -r12[1]/d12_c -r13[1]/d13_c;
  f[2] = -r12[2]/d12_c -r13[2]/d13_c;
  
  f[3] =  r12[0]/d12_c -r23[0]/d23_c;
  f[4] =  r12[1]/d12_c -r23[1]/d23_c;
  f[5] =  r12[2]/d12_c -r23[2]/d23_c;

  if( df != NULL ){
    //We are also requesting the Jacobian of the force!   
    //fifth power of distance is needed
    double d12_p = d12_c*d12*d12;
    double d13_p = d13_c*d13*d13;
    double d23_p = d23_c*d23*d23;

    //Use macros
    //eventually do subexpression elimination
    for( int i=0; i<3; i++ ){
      for( int j=0; j<3; j++){
        df[IDX(i  ,j  ,6,6)] =  HELPER_JACOBIAN(r12,d12_c,d12_p,i,j) + 2*HELPER_JACOBIAN(r13,d13_c,d13_p,i,j);
        df[IDX(i  ,j+3,6,6)] = -HELPER_JACOBIAN(r12,d12_c,d12_p,i,j) +   HELPER_JACOBIAN(r13,d13_c,d13_p,i,j);
        
        df[IDX(i+3,j  ,6,6)] = -HELPER_JACOBIAN(r12,d12_c,d12_p,i,j) +   HELPER_JACOBIAN(r23,d23_c,d23_p,i,j);
        df[IDX(i+3,j+3,6,6)] =  HELPER_JACOBIAN(r12,d12_c,d12_p,i,j) + 2*HELPER_JACOBIAN(r23,d23_c,d23_p,i,j);
      }
    }
  }

  return 0; //everything went well!
}

double hamiltonian( double *x ){
  //Separation vectors
  double r12[3] = {   x[0] -   x[3],   x[1] -   x[4],   x[2] -   x[5] }; //r1 - r2
  double r13[3] = { 2*x[0] +   x[3], 2*x[1] +   x[4], 2*x[2] +   x[5] }; //r1 - r3
  double r23[3] = {   x[0] + 2*x[3],   x[1] + 2*x[4],   x[2] + 2*x[5] }; //r2 - r3

  //Compute magnitudes
  double d12 = VECNORM(r12);
  double d13 = VECNORM(r13);
  double d23 = VECNORM(r23);
 
  //Get pointers to momentum
  double *p1 = &x[6];
  double *p2 = &x[9];

  //Compute energy
  double E = DOT3D( p1, p1 ) + DOT3D( p1, p2 ) + DOT3D( p2, p2 ) - 1.0/d12 - 1.0/d13 - 1.0/d23;
  return E;
}

int rescale( double *x ){
  //rescale state such that H=-1. If H>0, return an error code
  double E = hamiltonian(x);
  if( E > 0 ){
    //positive energy. Not a bound state. 
    return -1;
  }

  double l = sqrt(fabs(E));
  double *q = &x[0];
  double *p = &x[6];

  for( int i=0; i<6; i++){
    q[i] = q[i]*l*l; //rescale q
    p[i] = p[i]/l;   //rescale p
  }
  x[12] = x[12]*l*l*l; //rescale T

  l = hamiltonian( x );
  // printf("energy went from %f to %f\n", E, l);

  return 0;
}
