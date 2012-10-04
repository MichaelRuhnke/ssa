// g2o - General Graph Optimization
// Copyright (C) 2011 M. Ruhnke, 
// 
// g2o is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// g2o is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "pcg_cuda.h"
// Utilities and system includes
#include <shrUtils.h>
#include <shrQATest.h>
#include <cutil_inline.h>
#include <cusparse.h>
#include <cublas.h>

namespace g2o {

  SolverCUDACG::SolverCUDACG() { };
  SolverCUDACG::~SolverCUDACG() { };

    bool SolverCUDACG::solve(double* A, int* ptr, int* indices, int cols, int nz, double* x, double* b)
    {
      cusparseHandle_t   handle;
      cusparseStatus_t   status;
      cusparseMatDescr_t descr;

      /* Get handle to the CUSPARSE context */
      status = cusparseCreate(&handle);
      /* Description of the A matrix*/
      status = cusparseCreateMatDescr(&descr); 
      cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_SYMMETRIC);
      cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);
      cusparseSetMatFillMode (descr, CUSPARSE_FILL_MODE_UPPER );

      if(status != 0)
        return false;
      /** copy ccs to gpu */
      //Sparse matrix
      double* cu_A_ccs;
      int *cu_col_ptr_ccs, *cu_row_idx_ccs;

      double *cu_x, *cu_b; 
      double *cu_p, *cu_Ap;
      double *cu_r;

      //alloc memory
      cutilSafeCall( cudaMalloc((void**)&cu_col_ptr_ccs, sizeof(int)*(cols+1)) );
      cutilSafeCall( cudaMalloc((void**)&cu_row_idx_ccs, sizeof(int)*nz) );
      cutilSafeCall( cudaMalloc((void**)&cu_A_ccs, sizeof(double)*nz) );
      cutilSafeCall( cudaMalloc((void**)&cu_x, sizeof(double)*cols) );  
      cutilSafeCall( cudaMalloc((void**)&cu_b, sizeof(double)*cols) );
      cutilSafeCall( cudaMalloc((void**)&cu_r, sizeof(double)*cols) );
      cutilSafeCall( cudaMalloc((void**)&cu_p, sizeof(double)*cols) );
      cutilSafeCall( cudaMalloc((void**)&cu_Ap, sizeof(double)*cols) );

      cudaMemcpy(cu_A_ccs, A, sizeof(double)*nz, cudaMemcpyHostToDevice);
      cudaMemcpy(cu_row_idx_ccs, indices, sizeof(int)*nz, cudaMemcpyHostToDevice);
      cudaMemcpy(cu_col_ptr_ccs, ptr, sizeof(int)*(cols+1), cudaMemcpyHostToDevice);
      cudaMemcpy(cu_x, x, sizeof(double)*cols, cudaMemcpyHostToDevice);
      cudaMemcpy(cu_b, b, sizeof(double)*cols, cudaMemcpyHostToDevice);

      /** Solve CG problem */
      int max_iter = cols;
      int iteration = 0;
      double tol = 1e-6;
      double alpha, beta;

      /**  r=b-A*x; */
        /** r = b; */
        cublasDcopy (cols, cu_b, 1, cu_r, 1); 
        /** r = -1.0 * A * x + 1.0 * r; (r = b previous step) */
        cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, cols, cols, -1.0, 
                          descr, cu_A_ccs, cu_col_ptr_ccs, cu_row_idx_ccs, cu_x, 1.0, cu_r);
        /** p = r; */
        cublasDcopy (cols, cu_r, 1, cu_p, 1); 
      /** rs_old = r'*r; */
      double rs_old = cublasDdot(cols, cu_r, 1, cu_r, 1);
      double rs_new = 0;

      while (rs_old > tol*tol && iteration <= max_iter) {
        iteration++;

        cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, cols, cols, 1.0, 
                        descr, cu_A_ccs, cu_col_ptr_ccs, cu_row_idx_ccs, cu_p, 0.0, cu_Ap);   /** Ap = 1.0 * A * p + 0.0 * y; */
          alpha = rs_old/cublasDdot (cols, cu_Ap, 1, cu_p, 1);  /** alpha=rs_old/(p'*Ap); */
          cublasDaxpy (cols, alpha, cu_p, 1, cu_x, 1);          /** x=x+alpha*p; */
          cublasDaxpy (cols, -alpha, cu_Ap, 1, cu_r, 1);        /** r=r-alpha*Ap; */
          rs_new = cublasDdot(cols, cu_r, 1, cu_r, 1);          /** rsnew=r'*r; */
        /** p=r+rsnew/rsold*p; */
          beta = rs_new/rs_old; /** beta = rsnew/rsold */
          cublasDscal(cols, beta, cu_p, 1); /** p=beta*p; */
          cublasDaxpy (cols, 1.0, cu_r, 1, cu_p, 1) ;/** p=r+p; */
        rs_old=rs_new;
      }
      shrLog("  iteration = %3d, residual = %e \n", iteration, sqrt(rs_old));
      cudaMemcpy(x, cu_x, sizeof(double)*cols, cudaMemcpyDeviceToHost);

      cudaFree(cu_col_ptr_ccs);
      cudaFree(cu_row_idx_ccs);
      cudaFree(cu_A_ccs);
      cudaFree(cu_x);
      cudaFree(cu_r);
      cudaFree(cu_p);
      cudaFree(cu_Ap);
      cudaFree(cu_b);
      return true;
    }

   bool SolverCUDACG::solve(float* A, int* ptr, int* indices, int cols, int nz, float* x, float* b)
    {
      cusparseHandle_t   handle;
      cusparseStatus_t   status;
      cusparseMatDescr_t descr;

      /* Get handle to the CUSPARSE context */
      status = cusparseCreate(&handle);
      /* Description of the A matrix*/
      status = cusparseCreateMatDescr(&descr); 
      cusparseSetMatType(descr,CUSPARSE_MATRIX_TYPE_SYMMETRIC);
      cusparseSetMatIndexBase(descr,CUSPARSE_INDEX_BASE_ZERO);
      cusparseSetMatFillMode (descr, CUSPARSE_FILL_MODE_UPPER );

      if(status != 0)
        return false;
      /** copy ccs to gpu */
      //Sparse matrix
      float* cu_A_ccs;
      int *cu_col_ptr_ccs, *cu_row_idx_ccs;

      float *cu_x, *cu_b; 
      float *cu_p, *cu_Ap;
      float *cu_r;

      //alloc memory
      cutilSafeCall( cudaMalloc((void**)&cu_col_ptr_ccs, sizeof(int)*(cols+1)) );
      cutilSafeCall( cudaMalloc((void**)&cu_row_idx_ccs, sizeof(int)*nz) );
      cutilSafeCall( cudaMalloc((void**)&cu_A_ccs, sizeof(float)*nz) );
      cutilSafeCall( cudaMalloc((void**)&cu_r, sizeof(float)*cols) );
      cutilSafeCall( cudaMalloc((void**)&cu_x, sizeof(float)*cols) );  
      cutilSafeCall( cudaMalloc((void**)&cu_b, sizeof(float)*cols) );
      cutilSafeCall( cudaMalloc((void**)&cu_p, sizeof(float)*cols) );
      cutilSafeCall( cudaMalloc((void**)&cu_Ap, sizeof(float)*cols) );

      cudaMemcpy(cu_A_ccs, A, sizeof(float)*nz, cudaMemcpyHostToDevice);
      cudaMemcpy(cu_row_idx_ccs, indices, sizeof(int)*nz, cudaMemcpyHostToDevice);
      cudaMemcpy(cu_col_ptr_ccs, ptr, sizeof(int)*(cols+1), cudaMemcpyHostToDevice);
      cudaMemcpy(cu_x, x, sizeof(float)*cols, cudaMemcpyHostToDevice);
      cudaMemcpy(cu_b, b, sizeof(float)*cols, cudaMemcpyHostToDevice);

      /** Solve CG problem */
      int max_iter = cols;
      int iteration = 0;
      float tol = 1e-8;
      float alpha, beta;

      /**  r=b-A*x; */
        /** r = b; */
        cublasScopy (cols, cu_b, 1, cu_r, 1); 
        /** r = -1.0 * A * x + 1.0 * r; (r = b previous step) */
        cusparseScsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, cols, cols, -1.0, 
                          descr, cu_A_ccs, cu_col_ptr_ccs, cu_row_idx_ccs, cu_x, 1.0, cu_r);
        /** p = r; */
        cublasScopy (cols, cu_r, 1, cu_p, 1); 
      /** rs_old = r'*r; */
      float rs_old = cublasSdot(cols, cu_r, 1, cu_r, 1);
      float rs_new = 0;

      while (rs_old > tol*tol && iteration <= max_iter) {
        iteration++;

        cusparseScsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, cols, cols, 1.0, 
                        descr, cu_A_ccs, cu_col_ptr_ccs, cu_row_idx_ccs, cu_p, 0.0, cu_Ap);   /** Ap = 1.0 * A * p + 0.0 * y; */
          alpha = rs_old/cublasSdot (cols, cu_Ap, 1, cu_p, 1);  /** alpha=rs_old/(p'*Ap); */
          cublasSaxpy (cols, alpha, cu_p, 1, cu_x, 1);          /** x=x+alpha*p; */
          cublasSaxpy (cols, -alpha, cu_Ap, 1, cu_r, 1);        /** r=r-alpha*Ap; */
          rs_new = cublasSdot(cols, cu_r, 1, cu_r, 1);          /** rsnew=r'*r; */
        /** p=r+rsnew/rsold*p; */
          beta = rs_new/rs_old; /** beta = rsnew/rsold */
          cublasSscal(cols, beta, cu_p, 1); /** p=beta*p; */
          cublasSaxpy (cols, 1.0, cu_r, 1, cu_p, 1) ;/** p=r+p; */
        rs_old=rs_new;
      }
      shrLog("  iteration = %3d, residual = %e \n", iteration, sqrt(rs_old));
      cudaMemcpy(x, cu_x, sizeof(float)*cols, cudaMemcpyDeviceToHost);

      cudaFree(cu_col_ptr_ccs);
      cudaFree(cu_row_idx_ccs);
      cudaFree(cu_A_ccs);
      cudaFree(cu_x);
      cudaFree(cu_r);
      cudaFree(cu_p);
      cudaFree(cu_Ap);
      cudaFree(cu_b);
      return true;
    }

}// end namespace


