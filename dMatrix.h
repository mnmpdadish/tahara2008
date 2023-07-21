#pragma once

#include "utilities.h"

// reminder:
// int * ptr;
// int const * ptrToConst;
// int * const constPtr;
// int const * const constPtrToConst;
// At first, I added const AFTER each *, but in the end, it is just more confusing than nothing.

unsigned int dgemm_(char const*, char const*, unsigned int const*, unsigned int const*, unsigned int const*,
                  double const*, double const*, unsigned int const*, double const*, 
                  unsigned int const*, double const*, double *, unsigned int const*);
unsigned int dswap_(unsigned int const*, double*, unsigned int const*, double*, unsigned int const*);
unsigned int dgetrf_(unsigned int const*, unsigned int const*, double const*, unsigned int const*,unsigned int*,unsigned int*);
unsigned int dgetri_(unsigned int const*, double*, unsigned int const*, unsigned int const*, double*, unsigned int const*, unsigned int*);
unsigned int dgemv_(char const*, unsigned int const*, unsigned int const*, double const*, double const*, unsigned int const*, double const*, unsigned int const*, double const*, double*, unsigned int const*);
double ddot_(unsigned int const*, double const*, unsigned int const*, double const*, unsigned int const*);
unsigned int dger_(unsigned int const*, unsigned int const*, double const*, double const*, unsigned int const*, double const*, unsigned int const*, double *, unsigned int const*);
unsigned int dcopy_(unsigned int const*, double const*, unsigned int const* , double*, unsigned int const*);
unsigned int daxpy_(unsigned int const*, double const*, double const*, unsigned int const*, double*, unsigned int const*);
unsigned int dscal_(unsigned int const*, double const*, double*, unsigned int const*);

unsigned int dskpfa_(char const*, char const*, unsigned int const*, double*, unsigned int const*, double*, unsigned int*, double*, unsigned int const*, unsigned int*);




typedef struct {
  unsigned int N;
  unsigned int capacity;
  double *data;
} dMatrix;

typedef dMatrix dVector;

unsigned int resize_dVector(dVector * x, unsigned int N) {
  while(N > x->capacity) x->data = (double*) realloc(x->data, (x->capacity *= 2) * sizeof(double));
  x->N=N;
  return 0;
}

unsigned int resize_dMatrix(dMatrix * A, unsigned int N) {
  while(N*N > A->capacity) {
    A->data = (double*) realloc(A->data, (A->capacity *= 2) * sizeof(double));
    //printf("new matrix size = %d \n",A->capacity); fflush(stdout);
  }
  A->N=N;
  return 0;
}

unsigned int reset_dVector(dVector * X) {
  unsigned int i;
  for(i=0; i<(X->N); i++) X->data[i]=0;
  return 0;
}

unsigned int reset_dMatrix(dMatrix * A) {
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) A->data[i]=0;
  return 0;
}

unsigned int init_dVector(dVector * x, unsigned int N) {
  x->capacity = INIT_CAPACITY;
  x->data = (double *) malloc(x->capacity * sizeof(double));
  x->N=N;
  resize_dVector(x,N);
  return 0;
}

unsigned int init_dMatrix(dMatrix * A, unsigned int N) {
  A->capacity = INIT_CAPACITY;
  A->data = (double *) malloc(A->capacity * sizeof(double));
  A->N=N;
  resize_dMatrix(A,N);
  return 0;
}

unsigned int free_dVector(dVector * x) {
  free(x->data);
  return 0;
}

unsigned int free_dMatrix(dMatrix * A) {
  //Array_double_free(&A->buffer);
  free(A->data);
  return 0;
}



//Y==X
unsigned int areEqual_dVector(dVector const * X, dVector const * Y) {
  if(X->N != Y->N) return 0;
  unsigned int i;
  for(i=0; i<X->N; i++) if(!doubleEqual(X->data[i], Y->data[i])) return 0;
  return 1;
}

//B==A
unsigned int areEqual_dMatrix(dMatrix const * A, dMatrix * B) {
  if(A->N != B->N) return 0;
  unsigned int i;
  for(i=0; i<(A->N*A->N); i++) if(!doubleEqual(A->data[i], B->data[i])) return 0;
  return 1;
}


//Y=X
unsigned int copy_dVector(dVector const * X, dVector * Y) {
  //Y->N=X->N;
  resize_dVector(Y,X->N);
  memcpy(Y->data, X->data, X->N * sizeof(double));
  return 0;
}

//B=A
unsigned int copy_dMatrix(dMatrix const * A, dMatrix * B) {
  //B->N=A->N;
  resize_dMatrix(B,A->N);
  memcpy(B->data, A->data, A->N*A->N * sizeof(double));
  return 0;
}


#define ELEM(mtx, i, j) (mtx->data[j * mtx->N + i])
#define ELEM_VAL(mtx, i, j) (mtx.data[j * mtx.N + i])
//B=A
//where dim(B)=NxN and dim(A)=(A->N)x(A->N)
unsigned int copySub_dMatrix(dMatrix const * A, dMatrix * B, unsigned int N) {
  assert(N>=0);
  unsigned int copyLength = (A->N < N) ? A->N : N; //choose the smallest dimension between A->N and N. 
  resize_dMatrix(B,N);
  
  unsigned int i=0;	
  for(i = 0; i < copyLength; i++) {
    memcpy(&ELEM(B,0,i), &ELEM(A,0,i), sizeof(double)*copyLength);
  }
  return 0;
}


//C=A*B
unsigned int dMatrixMatrixMultiplication(dMatrix const * A, dMatrix const * B, dMatrix * C) {
  assert(A->N == B->N);
  resize_dMatrix(C,A->N);
  unsigned int N=A->N;
  double one=1.0;
  double zero=0.0;
  char no = 'n';
  dgemm_(&no,&no,&N,&N,&N, &one, A->data, &N, B->data, &N, &zero, C->data, &N); 
  return 0;
}

unsigned int dMatrixSwapCols(dMatrix * A, unsigned int row1, unsigned int row2) {
  unsigned int N=A->N;
  assert(row1 < N);
  assert(row2 < N);
  unsigned int one=1;
  dswap_(&N, &A->data[row1*A->N], &one, &A->data[row2*A->N], &one);
  return 0;
}

unsigned int dMatrixSwapRows(dMatrix * A, unsigned int col1, unsigned int col2) {
  unsigned int N=A->N;
  assert(col1 < N);
  assert(col2 < N);
  dswap_(&N, &A->data[col1], &N, &A->data[col2], &N);
  return 0;
}

void swap_double(double * a, double * b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

void swap_ptrOfDoubles(double * a, double * b)
{
    double * temp = a;
    a = b;
    b = temp;
}


//A^T
unsigned int transpose_dMatrix(dMatrix * A) {
  unsigned int i,j;
  for(i=0;i<A->N;i++)
    for(j=0;j<i;j++){
      swap_double(&ELEM(A,i,j),&ELEM(A,j,i));
    }
  return 0;
}


unsigned int print_dMatrix(dMatrix const * A) {
  if (!A){ 
    printf("oups.\n");
    return -1;
  }
  
  unsigned int i,j;
  double val;//, tol=0.00001;
  for (i = 0; i < A->N; i++) {
    for (j = 0; j < A->N; j++) {
      val = ELEM(A, i, j);
      //if(fabs(val)>tol) 
      printf("% 4.4f ", val);
      //else printf("   .   ");
    }
    printf("\n");
  }
  fflush(stdout);
  return 0;
}


unsigned int print_dVector(dVector const * X) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  unsigned int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.3f ", X->data[i]);
  }
  printf("\n");

  fflush(stdout);
  return 0;
}


unsigned int print_dVectorFactor(dVector const * X, double f) {
  if (!X){ 
    printf("oups.\n");
    return -1;
  }
  unsigned int i;
  for (i = 0; i < X->N; i++) {
    printf("% 6.3f ", f*X->data[i]);
  }
  printf("\n");

  fflush(stdout);
  return 0;
}


// A=>A^-1
void invert_dMatrix(dMatrix * A) {
  unsigned int INFO1=0;
  unsigned int INFO2=0;
  //unsigned int IPIV[A->N];
  unsigned int nEntry=A->N*A->N;
  //double WORK[nEntry]; // WARNING, FOR TOO BIG MATRICES, this stack declaration causes segmentation fault (happenned above 1200*1200 matrices for my computer).
  unsigned int * IPIV;
  double * WORK;
  WORK = (double *) malloc(nEntry*sizeof(double));
  IPIV = (unsigned int *) malloc(A->N*sizeof(unsigned int));
  
  dgetrf_(&A->N, &A->N, &A->data[0], &A->N, &IPIV[0], &INFO1);
  dgetri_(&A->N, A->data, &A->N, &IPIV[0], &WORK[0], &nEntry, &INFO2);
  if( !(INFO1 == 0) || !(INFO2 == 0) ) {
    printf( "The algorithm failed to invert the dMatrix. %d %d\n", INFO1, INFO2);
    exit( 1 );
  }
  free(WORK);
  free(IPIV);
}

// Y=A*X
void dMatrixVectorProduct(dMatrix const*A, dVector const*X, double const factor, dVector *Y) {
  unsigned int N=A->N;
  assert(N == X->N);
  //assert(N == Y->N);
  resize_dVector(Y,N);
  //double one=1.0;
  double zero=0.0;
  char no = 'n';
  unsigned int inc=1;
  dgemv_(&no, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}

void dVectorMatrixProduct(dVector const*X, dMatrix const*A, double const factor, dVector *Y) {
  unsigned int N=A->N;
  assert(N == X->N);
  //assert(N == Y->N);
  resize_dVector(Y,N);
  //double one=1.0;
  double zero=0.0;
  char yes = 't';
  unsigned int inc=1;
  dgemv_(&yes, &N, &N, &factor, A->data, &N, X->data, &inc, &zero, Y->data, &inc); 
}


// return X.Y
double dScalarProduct(dVector const*X, dVector const*Y) {
  assert(X->N==Y->N);
  unsigned int inc=1;
  return ddot_(&X->N, X->data, &inc, Y->data, &inc);
}

void scale_dVector(dVector *X, double const scal){
  unsigned int inc=1;
  dscal_(&X->N, &scal, X->data, &inc);
}



// Suppose a matrix A composed of the block matrices Aij:
// A = [ A11 A12 ]
//     [ A21 A22 ]
//
// which is the inverse of D, such that:
// [ D11 D12 ]  [ A11 A12 ] = [ 1 0 ]
// [ D21 D22 ]  [ A11 A12 ]   [ 0 1 ]
//
// then we can prove directly that the Schur complement is:
// D11^(-1) = A11 - A12 A22^(-1) A21
//
// At the output of the function S = D11^-1 (the Schur complement).
// Here, dimensions are:
// dim(A) = NxN, 
// dim(A11)=(N-1)x(N-1), dim(A12)=1x(N-1), dim(A21)=(N-1)xN, dim(A22)=1x1
// Same for Dij and dim(S)=(N-1)x(N-1).
unsigned int dSchurComplement(dMatrix const*A, dMatrix *S) {
  //printf("before dSchur\n"); fflush(stdout);
  unsigned int inc = 1;
  unsigned int N = A->N;
  unsigned int Nm1 = N-1;
  //S->N=Nm1;
  //assert(S->N==N-1);
  //assert(Nm1>=0);
  if(Nm1==0) resize_dMatrix(S,0);
  else{
    copySub_dMatrix(A,S,Nm1);

    double factor = -1./ELEM(A,Nm1,Nm1);
    // this next line does S = S - A12 A22^(-1) A21;
    dger_(&S->N, &S->N, &factor, &ELEM(A,0,Nm1), &inc, &ELEM(A,Nm1,0), &A->N, S->data, &S->N);
  }
  //printf("after dSchur\n"); fflush(stdout);
  return 0;
}



unsigned int dAddRowColToInverse(dMatrix const*A, dVector const*Rtilde, dVector const*Qtilde, double const sTilde, dMatrix *Ap1) {
  //printf("before dAddRow\n"); fflush(stdout);
  unsigned int inc = 1;
  unsigned int N = A->N;
  assert(N==Rtilde->N);
  assert(N==Qtilde->N);
  unsigned int Np1 = N+1;
  copySub_dMatrix(A,Ap1,Np1);

  double factor = 1.0/sTilde;
  // got pAccept now.
  dger_(&N, &N, &factor, Qtilde->data, &inc, Rtilde->data, &inc, Ap1->data, &Np1);
  
  dcopy_(&N, Rtilde->data, &inc, &ELEM(Ap1,N,0), &Np1);
  dcopy_(&N, Qtilde->data, &inc, &ELEM(Ap1,0,N), &inc);
  ELEM(Ap1,N,N)=sTilde;
  //printf("after dAddRow\n"); fflush(stdout);
  return 0;
}



unsigned int dCopyColIntoVector(dMatrix const*A, dVector *Col, unsigned int index) {
  unsigned int N = A->N;
  resize_dVector(Col,N);
  //assert(N==Row->N);
  memcpy(Col->data, &ELEM(A,0,index), sizeof(double)*N);
  return 0;
}

unsigned int dCopyRowIntoVector(dMatrix const*A, dVector *Row, unsigned int index) {
  unsigned int inc = 1;
  unsigned int N = A->N;
  resize_dVector(Row,N);
  //assert(N==Row->N);
  dcopy_(&N, &ELEM(A,index,0), &N, Row->data, &inc);
  return 0;
}

unsigned int dAddOneElementToInverse(dMatrix *A, dVector const*Row, dVector const*Col, double const value) {
  unsigned int inc = 1;
  unsigned int N = A->N;
//  assert(N==Row->N);
//  assert(N==Col->N);
  dger_(&N, &N, &value, Col->data, &inc, Row->data, &inc, A->data, &N);
  
  return 0;
}


double dPfaffian(dMatrix *A) {
  
  // Ref Tahara, Imada 2008, B.11
  
  char up = 'U';
  char method = 'P';
  
  unsigned int N = A->N;
  double pfaffian;
  
  unsigned int LWORK = N*N;
  unsigned int INFO = 0;
  
  double * WORK;
  unsigned int * IWORK;
  
  dMatrix B;
  init_dMatrix(&B,0);
  copy_dMatrix(A,&B);
  
  
  WORK  = (double *) malloc(LWORK*sizeof(double));
  IWORK = (unsigned int *) malloc(N*sizeof(unsigned int));
  
  dskpfa_(&up, &method, &N, &B.data[0], &N, &pfaffian,  &IWORK[0], &WORK[0], &LWORK, &INFO);
  if( !(INFO == 0) ) {
    printf( "The algorithm failed to calculate the Pfaffian. %d\n", INFO);
    exit( 1 );
  }
    
  free(WORK);
  free(IWORK);
  return pfaffian;
}



unsigned int dUpdateHopping(double pf_A, dMatrix const *A_inv, dMatrix *B_inv, dMatrix const*bm, int alpha) {

  // Ref Tahara, Imada 2008, B.14
  
  copy_dMatrix(A_inv,B_inv);
  unsigned int N = B_inv->N;
  
  dVector X;
  init_dVector(&X,N);  
  dMatrixVectorProduct(B_inv,bm, 1.0,&X);
  
/*  implementation 1:
  int j=0, i=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      if((i!=alpha) && (j!=alpha)) {
        ELEM(B_inv, i, j) += (1./(X.data[alpha])) * ( -X.data[i]*ELEM(A_inv, alpha, j) + X.data[j]*ELEM(A_inv, alpha, i) ); 
      }
      else{
        ELEM(B_inv, i, j) *= (1./(X.data[alpha]));       
      }
    }
  }
*/
  
/*  implementation 2:
  int j=0, i=0;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      //if((i!=alpha) && (j!=alpha)) {
        ELEM(B_inv, i, j) += (1./(X.data[alpha])) * ( -X.data[i]*ELEM(A_inv, alpha, j) + X.data[j]*ELEM(A_inv, alpha, i) ); 

      //kronecker i,alpha:
      if(i==alpha) {
        ELEM(B_inv, i, j) += (1./(X.data[alpha])) * ( ELEM(A_inv, i, j) ); 
      }
      //kronecker j,alpha:
      if(j==alpha) {
        ELEM(B_inv, i, j) += (1./(X.data[alpha])) * ( ELEM(A_inv, i, j) ); 
      }
    }
  }
*/
  
  //implementation 3, no loop:
  unsigned int inc = 1;
  double factor_pos = (1./(X.data[alpha]));
  double factor_neg =-(1./(X.data[alpha]));
  
  // B^-1 = B^-1 - 1./(X[alpha])) * X*Y^T             where Y is the m column of A_inv
  dger_(&N, &N, &factor_neg, X.data, &inc, &A_inv->data[alpha], &N, B_inv->data, &N);
  
  // B^-1 = B^-1 + 1./(X[alpha])) * Y*X^T
  dger_(&N, &N, &factor_pos, &A_inv->data[alpha], &N, X.data, &inc, B_inv->data, &N);
  
  // adding column alpha of A_inv to B_inv
  daxpy_(&N, &factor_pos, &A_inv->data[alpha], &N, &B_inv->data[alpha], &N);
  
  // adding row alpha of A_inv to B_inv
  daxpy_(&N, &factor_pos, &A_inv->data[N*alpha], &inc, &B_inv->data[N*alpha], &inc);

  //if(verbose) {printf("\nX.data[alpha]= % 4.5f\n", X.data[alpha]);}
      
  return 0;
}


