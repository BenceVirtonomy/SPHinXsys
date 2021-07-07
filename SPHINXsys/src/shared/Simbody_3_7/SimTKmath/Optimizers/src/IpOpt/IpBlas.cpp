// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpBlas.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpBlas.hpp"

#define WASM_BUILD
#ifndef WASM_BUILD

#if SimTK_DEFAULT_PRECISION==1 // float
#define DCOPY   scopy_
#define DSCAL   sscal_
#define DAXPY   saxpy_
#define DDOT    sdot_
#define DNRM2   snrm2_
#define DASUM   sasum_
#define DGEMV   sgemv_
#define DSYMV   ssymv_
#define DGEMM   sgemm_
#define DSYRK   ssyrk_
#define DTRSM   strsm_
#define IDAMAX  isamax_
#else // double
#define DCOPY   dcopy_
#define DSCAL   dscal_
#define DAXPY   daxpy_
#define DDOT    ddot_
#define DNRM2   dnrm2_
#define DASUM   dasum_
#define DGEMV   dgemv_
#define DSYMV   dsymv_
#define DGEMM   dgemm_
#define DSYRK   dsyrk_
#define DTRSM   dtrsm_
#define IDAMAX  idamax_
#endif

// Prototypes for the BLAS routines
extern "C"
{
  /** BLAS Fortran function DDOT */
  double ddot_(ipfint *n, const double *x, ipfint *incX,
                             const double *y, ipfint *incY);
  /** BLAS Fortran function DNRM2 */
  double dnrm2_(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran function DASUM */
  double dasum_(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran function IDAMAX */
  ipfint idamax_(ipfint *n, const double *x, ipfint *incX);
  /** BLAS Fortran subroutine DCOPY */
  void dcopy_(ipfint *n, const double *x, ipfint *incX,
                             double *y, ipfint *incY);
  /** BLAS Fortran subroutine DAXPY */
  void daxpy_(ipfint *n, const double *alpha, const double *x,
                             ipfint *incX, double *y, ipfint *incY);
  /** BLAS Fortran subroutine DSCAL */
  void dscal_(ipfint *n, const double *alpha, const double *x,
                             ipfint *incX);

  void dgemv_(char* trans, ipfint *m, ipfint *n,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *x, ipfint *incX, const double *beta,
                             double *y, ipfint *incY, int trans_len);
  void dsymv_(char* uplo, ipfint *n,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *x, ipfint *incX, const double *beta,
                             double *y, ipfint *incY, int uplo_len);
  void dgemm_(char* transa, char* transb,
                             ipfint *m, ipfint *n, ipfint *k,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *b, ipfint *ldb, const double *beta,
                             double *c, ipfint *ldc,
                             int transa_len, int transb_len);
  void dsyrk_(char* uplo, char* trans, ipfint *n, ipfint *k,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *beta, double *c, ipfint *ldc,
                             int uplo_len, int trans_len);
  void dtrsm_(char* side, char* uplo, char* transa, char* diag,
                             ipfint *m, ipfint *n,
                             const double *alpha, const double *a, ipfint *lda,
                             const double *b, ipfint *ldb,
                             int side_len, int uplo_len,
                             int transa_len, int diag_len);

    /** BLAS Fortran function DDOT */
  float sdot_(ipfint *n, const float *x, ipfint *incX,
                             const float *y, ipfint *incY);
  /** BLAS Fortran function DNRM2 */
  float snrm2_(ipfint *n, const float *x, ipfint *incX);
  /** BLAS Fortran function DASUM */
  float sasum_(ipfint *n, const float *x, ipfint *incX);
  /** BLAS Fortran function IDAMAX */
  ipfint isamax_(ipfint *n, const float *x, ipfint *incX);
  /** BLAS Fortran subroutine DCOPY */
  void scopy_(ipfint *n, const float *x, ipfint *incX,
                             float *y, ipfint *incY);
  /** BLAS Fortran subroutine DAXPY */
  void saxpy_(ipfint *n, const float *alpha, const float *x,
                             ipfint *incX, float *y, ipfint *incY);
  /** BLAS Fortran subroutine DSCAL */
  void sscal_(ipfint *n, const float *alpha, const float *x,
                             ipfint *incX);

  void sgemv_(char* trans, ipfint *m, ipfint *n,
                             const float *alpha, const float *a, ipfint *lda,
                             const float *x, ipfint *incX, const float *beta,
                             float *y, ipfint *incY, int trans_len);
  void ssymv_(char* uplo, ipfint *n,
                             const float *alpha, const float *a, ipfint *lda,
                             const float *x, ipfint *incX, const float *beta,
                             float *y, ipfint *incY, int uplo_len);
  void sgemm_(char* transa, char* transb,
                             ipfint *m, ipfint *n, ipfint *k,
                             const float *alpha, const float *a, ipfint *lda,
                             const float *b, ipfint *ldb, const float *beta,
                             float *c, ipfint *ldc,
                             int transa_len, int transb_len);
  void ssyrk_(char* uplo, char* trans, ipfint *n, ipfint *k,
                             const float *alpha, const float *a, ipfint *lda,
                             const float *beta, float *c, ipfint *ldc,
                             int uplo_len, int trans_len);
  void strsm_(char* side, char* uplo, char* transa, char* diag,
                             ipfint *m, ipfint *n,
                             const float *alpha, const float *a, ipfint *lda,
                             const float *b, ipfint *ldb,
                             int side_len, int uplo_len,
                             int transa_len, int diag_len);
}

#endif // WASM_BUILD

namespace SimTKIpopt
{
#ifndef HAVE_CBLAS
  /* Interface to FORTRAN routine DDOT. */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDdot called"));
    return 0;
#elif WASM_BUILD
    ipfint n=size, INCX=incX, INCY=incY;

    return DDOT(&n, x, &INCX, y, &INCY);
#endif
  }

  /* Interface to FORTRAN routine DNRM2. */
  Number IpBlasDnrm2(Index size, const Number *x, Index incX)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDnrm2 called"));
    return 0;
#elif WASM_BUILD
    ipfint n=size, INCX=incX;

    return DNRM2(&n, x, &INCX);
#endif
  }

  /* Interface to FORTRAN routine DASUM. */
  Number IpBlasDasum(Index size, const Number *x, Index incX)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDasum called"));
    return 0;
#elif WASM_BUILD
    ipfint n=size, INCX=incX;

    return DASUM(&n, x, &INCX);
#endif
  }

  /* Interface to FORTRAN routine DASUM. */
  Index IpBlasIdamax(Index size, const Number *x, Index incX)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasIdamax called"));
    return 0;
#elif WASM_BUILD
    ipfint n=size, INCX=incX;

    return (Index) IDAMAX(&n, x, &INCX);
#endif
  }

  /* Interface to FORTRAN routine DCOPY. */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y, Index incY)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDcopy called"));
#elif WASM_BUILD
    ipfint N=size, INCX=incX, INCY=incY;

    DCOPY(&N, x, &INCX, y, &INCY);
#endif
  }

  /* Interface to FORTRAN routine DAXPY. */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX, Number *y,
                   Index incY)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDaxpy called"));
#elif WASM_BUILD
    ipfint N=size, INCX=incX, INCY=incY;

    DAXPY(&N, &alpha, x, &INCX, y, &INCY);
#endif
  }

  /* Interface to FORTRAN routine DSCAL. */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDscal called"));
#elif WASM_BUILD
    ipfint N=size, INCX=incX;

    DSCAL(&N, &alpha, x, &INCX);
#endif
  }

  void IpBlasDgemv(bool trans, Index nRows, Index nCols, Number alpha,
                   const Number* A, Index ldA, const Number* x,
                   Index incX, Number beta, Number* y, Index incY)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDgemv called"));
#elif WASM_BUILD
    ipfint M=nCols, N=nRows, LDA=ldA, INCX=incX, INCY=incY;

    char TRANS;
    if (trans) {
      TRANS = 'T';
    }
    else {
      TRANS = 'N';
    }

    DGEMV(&TRANS, &M, &N, &alpha, A, &LDA, x,
          &INCX, &beta, y, &INCY, 1);
#endif
  }

  void IpBlasDsymv(Index n, Number alpha, const Number* A, Index ldA,
                   const Number* x, Index incX, Number beta, Number* y,
                   Index incY)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDsymv called"));
#elif WASM_BUILD
    ipfint N=n, LDA=ldA, INCX=incX, INCY=incY;

    char UPLO='L';

    DSYMV(&UPLO, &N, &alpha, A, &LDA, x,
          &INCX, &beta, y, &INCY, 1);
#endif
  }

  void IpBlasDgemm(bool transa, bool transb, Index m, Index n, Index k,
                   Number alpha, const Number* A, Index ldA, const Number* B,
                   Index ldB, Number beta, Number* C, Index ldC)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDgemm called"));
#elif WASM_BUILD
    ipfint M=m, N=n, K=k, LDA=ldA, LDB=ldB, LDC=ldC;

    char TRANSA;
    if (transa) {
      TRANSA = 'T';
    }
    else {
      TRANSA = 'N';
    }
    char TRANSB;
    if (transb) {
      TRANSB = 'T';
    }
    else {
      TRANSB = 'N';
    }

    DGEMM(&TRANSA, &TRANSB, &M, &N, &K, &alpha, A, &LDA,
          B, &LDB, &beta, C, &LDC, 1, 1);
#endif
  }

  void IpBlasDsyrk(bool trans, Index ndim, Index nrank,
                   Number alpha, const Number* A, Index ldA,
                   Number beta, Number* C, Index ldC)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDsyrk called"));
#elif WASM_BUILD
    ipfint N=ndim, K=nrank, LDA=ldA, LDC=ldC;

    char UPLO='L';
    char TRANS;
    if (trans) {
      TRANS = 'T';
    }
    else {
      TRANS = 'N';
    }

    DSYRK(&UPLO, &TRANS, &N, &K, &alpha, A, &LDA,
          &beta, C, &LDC, 1, 1);
#endif
  }

  void IpBlasDtrsm(bool trans, Index ndim, Index nrhs, Number alpha,
                   const Number* A, Index ldA, Number* B, Index ldB)
  {
#ifdef WASM_BUILD
    throw std::runtime_error(std::string("SimTKIpopt::IpBlasDtrsm called"));
#elif WASM_BUILD
    ipfint M=ndim, N=nrhs, LDA=ldA, LDB=ldB;

    char SIDE = 'L';
    char UPLO = 'L';
    char TRANSA;
    if (trans) {
      TRANSA = 'T';
    }
    else {
      TRANSA = 'N';
    }
    char DIAG = 'N';

    DTRSM(&SIDE, &UPLO, &TRANSA, &DIAG, &M, &N,
          &alpha, A, &LDA, B, &LDB, 1, 1, 1, 1);
#endif
  }

#else
  /* Interface to CBLAS routine DDOT. */
  Number IpBlasDdot(Index size, const Number *x, Index incX, const Number *y,
                    Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DNRM2. */
  Number IpBlasDnrm2(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DASUM. */
  Number IpBlasDasum(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DASUM. */
  Index IpBlasIdamax(Index size, const Number *x, Index incX)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DCOPY. */
  void IpBlasDcopy(Index size, const Number *x, Index incX, Number *y, Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DAXPY. */
  void IpBlasDaxpy(Index size, Number alpha, const Number *x, Index incX, Number *y,
                   Index incY)
  {
    Not Implemented Yet!
  }

  /* Interface to CBLAS routine DSCAL. */
  void IpBlasDscal(Index size, Number alpha, Number *x, Index incX)
  {
    Not Implemented Yet!
  }

#endif

} // namespace Ipopt
