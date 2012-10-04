// g2o - General Graph Optimization
// Copyright (C) 2011 M. Ruhnke, H. Strasdat
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

#ifndef LINEAR_SOLVER_CUDA_GRADIENT_DESCENT_H
#define LINEAR_SOLVER_CUDA_GRADIENT_DESCENT_H

#include "g2o/core/linear_solver.h"
#include "g2o/core/batch_stats.h"

#include <vector>
#include <utility>
#include<Eigen/Core>
#include<Eigen/Cholesky>

// Utilities and system includes
#include <shrUtils.h>
#include <shrQATest.h>
#include <cutil_inline.h>
#include <cusparse.h>
#include <cublas.h>

#include "pcg_cuda.h"

#include "EXTERNAL/g2o/g2o/stuff/timeutil.h"

namespace g2o {
  using namespace std;
  /**
   * \brief linear solver using PCG, pre-conditioner is block Jacobi
   */
  template <typename MatrixType>
  class LinearSolverCUDACG : public LinearSolver<MatrixType>
  {
    public:
      LinearSolverCUDACG() :
      LinearSolver<MatrixType>(), _doublePrecision(true)
      {
        // This will pick the best possible CUDA capable device
        int argc = 0;
        char** argv = 0;
        cudaDeviceProp deviceProp;
        int devID = cutilChooseCudaDevice(argc, argv);
        if (devID < 0) {
            printf("exiting...\n");
            cutilExit(argc, argv);
            exit(EXIT_SUCCESS);
        }
        cutilSafeCall( cudaGetDeviceProperties(&deviceProp, devID) );
    
        // Statistics about the GPU device
        //printf("> GPU device has %d Multi-Processors, SM %d.%d compute capabilities\n\n", 
        //deviceProp.multiProcessorCount, deviceProp.major, deviceProp.minor);
    
        int version = (deviceProp.major * 0x10 + deviceProp.minor);
        if(version < 0x11) 
        {
          const char * sSDKname   = "conjugateGradientPrecond";
          printf("%s: requires a minimum CUDA compute 1.1 capability\n", sSDKname);
          cutilDeviceReset();
          shrQAFinishExit(argc, (const char **)argv, QA_PASSED);
        }

      }

      virtual ~LinearSolverCUDACG()
      {
      cutilDeviceReset();
       if(_colum_count > 0)
          delete _colum_ptr;
       if(_non_zero_entries > 0)
        delete _row_indices;
      }

      typedef std::vector< MatrixType, Eigen::aligned_allocator<MatrixType> > MatrixVector;
      typedef std::vector< const MatrixType* > MatrixPtrVector;

      virtual bool init()
      {
        _colum_count = 0;
        _non_zero_entries = 0;
        _colum_ptr = NULL;
        _row_indices = NULL;
        return true;
      }

      /* checkStatus: concise method for verifying CUDA return status */
      int checkStatus ( cusparseStatus_t status, char *msg )
      {
          if ( status != CUSPARSE_STATUS_SUCCESS ) {
              fprintf (stderr, "!!!! CUSPARSE %s ERROR \n", msg);
              cerr << status << endl;
              return 1;
          }
          return 0;
      }

      protected:
        int  _colum_count;
        int  _non_zero_entries;
        int* _colum_ptr;
        int* _row_indices;
        cusparseHandle_t   _handle;
        cusparseStatus_t   _status;
        cusparseMatDescr_t _descr;

        int _max_iter ;
        int _iteration;
        double _tol;
        bool _doublePrecision;

      public:

      void setDoublePrecision(bool value){
        _doublePrecision = value;
      }


      /* genICP: Generate the Incomplete Cholesky Preconditioner for a symmetric tridiagonal.  
        Follows description from Golub & Van Loan, "Matrix Computations 3rd Ed.", section 10.3.2 */
      void genICP ( int *rowPtrs, float *vals, int N, int *colIndsICP, int *rowPtrsICP, float *valsICP )
      {
          // Define a lower triangular banded matrix with 2 bands.
          rowPtrsICP[0] = 0;
          colIndsICP[0] = 0;
          int inz = 1;
          for ( int k=1; k<N; k++ ) {
              rowPtrsICP[k] = inz;
              for ( int j=k-1; j<=k; j++ ) {
                  colIndsICP[inz] = j;
                  inz++;
              }
          }
          rowPtrsICP[N] = inz;
      
          // copy A into H
          valsICP[0] = vals[0];
          for ( int k=1; k<N; k++ ) {
              valsICP[rowPtrsICP[k]] = vals[rowPtrs[k]];
              valsICP[rowPtrsICP[k]+1] = vals[rowPtrs[k]+1];    
          }
      
          // construct H
          for ( int k=0; k<N; k++ ) {
              valsICP[rowPtrsICP[k+1]-1] = sqrt(valsICP[rowPtrsICP[k+1]-1]);
              if ( k < N-1 ) {
                  valsICP[rowPtrsICP[k+1]] /= valsICP[rowPtrsICP[k+1]-1];      
                  valsICP[rowPtrsICP[k+1]+1] -= valsICP[rowPtrsICP[k+1]]*valsICP[rowPtrsICP[k+1]];
              }
          }
      
          return;
      }


      /* genICP: Generate the Incomplete Cholesky Preconditioner for a symmetric tridiagonal.  
        Follows description from Golub & Van Loan, "Matrix Computations 3rd Ed.", section 10.3.2 */
      void genICPD ( int *rowPtrs, double *vals, int N, int *colIndsICP, int *rowPtrsICP, double *valsICP )
      {
          // Define a lower triangular banded matrix with 2 bands.
          rowPtrsICP[0] = 0;
          colIndsICP[0] = 0;
          int inz = 1;
          for ( int k=1; k<N; k++ ) {
              rowPtrsICP[k] = inz;
              for ( int j=k-1; j<=k; j++ ) {
                  colIndsICP[inz] = j;
                  inz++;
              }
          }
          rowPtrsICP[N] = inz;
      
          // copy A into H
          valsICP[0] = vals[0];
          for ( int k=1; k<N; k++ ) {
              valsICP[rowPtrsICP[k]] = vals[rowPtrs[k]];
              valsICP[rowPtrsICP[k]+1] = vals[rowPtrs[k]+1];    
          }
      
          // construct H
//           for ( int k=0; k<N; k++ ) {
//               valsICP[rowPtrsICP[k+1]-1] = sqrt(valsICP[rowPtrsICP[k+1]-1]);
//               if ( k < N-1 ) {
//                   valsICP[rowPtrsICP[k+1]] /= valsICP[rowPtrsICP[k+1]-1];      
//                   valsICP[rowPtrsICP[k+1]+1] -= valsICP[rowPtrsICP[k+1]]*valsICP[rowPtrsICP[k+1]];
//               }
//           }
      
          return;
      }

      bool solve(const SparseBlockMatrix<MatrixType>& A, double* x, double* b)
      {
        if(_doublePrecision)
        {
          return solveDoublePrecision(A,x,b);
        } else {
          return solveSinglePrecision(A,x,b);
        }

        //return solveDoublePrecisionPrecond(A,x,b);
      }

        bool solveSinglePrecision(const SparseBlockMatrix<MatrixType>& A, double* x, double* b)
      {

        /**  initialize ccs colum pointers for matrix A */
        if(_colum_count < A.cols()){
          delete _colum_ptr;
          _colum_ptr = (int*)malloc(sizeof(int)*(A.cols()+1));
        }
        _colum_count = A.cols();

        /** initialize ccs row indices for matrix A */
        if(_non_zero_entries < (int) A.nonZeros()){
          delete _row_indices;
          _row_indices = (int*)malloc(sizeof(int)*(A.nonZeros()));
        }
        _non_zero_entries = A.nonZeros();

        // copy Data from A to A_float
        double *A_ccs = NULL;
        A_ccs = (double*)malloc(sizeof(double)*_non_zero_entries);                            // csr values for matrix A
        _non_zero_entries = A.fillCCS(_colum_ptr, _row_indices, A_ccs, false);

        //Convert to float (My nvidia card is too old for double precision...:( )
        float* A_ccs_float = (float*)malloc(sizeof(float)*_non_zero_entries);                            // csr values for matrix A
        for (int i = 0; i < _non_zero_entries; i++) {
          A_ccs_float[i] = (float) A_ccs[i];
        }

        /** initalize and convert to float (b and x) */
        float* x_float = (float*)malloc(sizeof(float)*_colum_count);
        float* b_float = (float*)malloc(sizeof(float)*_colum_count);
       //setting x and b up
        for (int i = 0; i < _colum_count; i++) {
          x_float[i] = 0.0;          // Initial approximation of solution
          b_float[i] = b[i];
        }

        SolverCUDACG cuda_cg;
        cuda_cg.solve(A_ccs_float, _colum_ptr, _row_indices, _colum_count ,_non_zero_entries, x_float, b_float);

        for (int i = 0; i < _colum_count; i++) {
          x[i] = (double) x_float[i];          // Initial approximation of solution
        }

       delete A_ccs;
       delete A_ccs_float;
       delete x_float;
       delete b_float;
       return true;
     }

        bool solveDoublePrecision(const SparseBlockMatrix<MatrixType>& A, double* x, double* b)
      {
       double timing = get_time();
        /**  initialize ccs colum pointers for matrix A */
        if(_colum_count < A.cols()){
          delete _colum_ptr;
          _colum_ptr = (int*)malloc(sizeof(int)*(A.cols()+1));
        }
        _colum_count = A.cols();

        /** initialize ccs row indices for matrix A */
        if(_non_zero_entries < (int) A.nonZeros()){
          delete _row_indices;
          _row_indices = (int*)malloc(sizeof(int)*(A.nonZeros()));
        }
        _non_zero_entries = A.nonZeros();
  
        // get data in CCS format
        double *A_ccs = NULL;
        A_ccs = (double*)malloc(sizeof(double)*_non_zero_entries);                            // csr values for matrix A
        _non_zero_entries = A.fillCCS(_colum_ptr, _row_indices, A_ccs, false);

        /** initalize x */
        for (int i = 0; i < _colum_count; i++) {
          x[i] = 0.0;
        }
        cerr << "mem alloc and CCS creation took " << (get_time() - timing)*1000 << "ms" << endl;

        SolverCUDACG cuda_cg;
        cuda_cg.solve(A_ccs, _colum_ptr, _row_indices, _colum_count ,_non_zero_entries, x, b);

       return true;
     }

        bool solveDoublePrecisionPrecond(const SparseBlockMatrix<MatrixType>& A, double* x, double* b)
      {
       double timing = get_time();
        /**  initialize ccs colum pointers for matrix A */
        if(_colum_count < A.cols()){
          delete _colum_ptr;
          _colum_ptr = (int*)malloc(sizeof(int)*(A.cols()+1));
        }
        _colum_count = A.cols();

        /** initialize ccs row indices for matrix A */
        if(_non_zero_entries < (int) A.nonZeros()){
          delete _row_indices;
          _row_indices = (int*)malloc(sizeof(int)*(A.nonZeros()));
        }
        _non_zero_entries = A.nonZeros();
  
        // get data in CCS format
        double *A_ccs = NULL;
        A_ccs = (double*)malloc(sizeof(double)*_non_zero_entries);                            // csr values for matrix A
        _non_zero_entries = A.fillCCS(_colum_ptr, _row_indices, A_ccs, false);

        /** initalize x */
        for (int i = 0; i < _colum_count; i++) {
          x[i] = 0.0;
        }
        cerr << "mem alloc and CCS creation took " << (get_time() - timing)*1000 << "ms" << endl;

        timing = get_time();
        /** copy ccs to gpu */
        // variable in cuda memory 
        int *cu_col_ptr_ccs, *cu_row_idx_ccs;
        double* cu_A_ccs;
        double *cu_x, *cu_b; 
        double *cu_p, *cu_Ap;
        double *cu_r;
        cutilSafeCall( cudaMalloc((void**)&cu_col_ptr_ccs, sizeof(int)*(_colum_count+1)) );
        cutilSafeCall( cudaMalloc((void**)&cu_row_idx_ccs, sizeof(int)*_non_zero_entries) );
        cutilSafeCall( cudaMalloc((void**)&cu_A_ccs, sizeof(double)*_non_zero_entries) );
        cutilSafeCall( cudaMalloc((void**)&cu_r, sizeof(double)*_colum_count) );
        cutilSafeCall( cudaMalloc((void**)&cu_x, sizeof(double)*_colum_count) );  
        cutilSafeCall( cudaMalloc((void**)&cu_b, sizeof(double)*_colum_count) );
        cutilSafeCall( cudaMalloc((void**)&cu_p, sizeof(double)*_colum_count) );
        cutilSafeCall( cudaMalloc((void**)&cu_Ap, sizeof(double)*_colum_count) );
        cerr << "GPU memory allocation took " << (get_time() - timing)*1000 << "ms" << endl;

        timing = get_time();
        cudaMemcpy(cu_col_ptr_ccs, _colum_ptr, sizeof(int)*(_colum_count+1), cudaMemcpyHostToDevice);
        cudaMemcpy(cu_row_idx_ccs, _row_indices, sizeof(int)*_non_zero_entries, cudaMemcpyHostToDevice);
        cudaMemcpy(cu_A_ccs, A_ccs, sizeof(double)*_non_zero_entries, cudaMemcpyHostToDevice);
        cudaMemcpy(cu_x, x, sizeof(double)*_colum_count, cudaMemcpyHostToDevice);
        cudaMemcpy(cu_b, b, sizeof(double)*_colum_count, cudaMemcpyHostToDevice);
        cerr << "GPU data transfer took " << (get_time() - timing)*1000 << "ms" << endl;

        timing = get_time();
      int *colIndsICP = NULL;
      int *rowPtrsICP = NULL;
      double *valsICP  = NULL;
  
      int nzICP = _non_zero_entries;
  
      valsICP = (double *) malloc (sizeof(double)*_non_zero_entries);
      colIndsICP = (int *) malloc (sizeof(int)*_non_zero_entries);
      rowPtrsICP = (int *) malloc (sizeof(int)*(_colum_count+1));
  
      // generate the Incomplete Cholesky factor H (lower triangular) for the matrix A.
      //genICPD ( _colum_ptr, A_ccs, _colum_count, colIndsICP, rowPtrsICP, valsICP );
      cerr << "Incomplete Cholesky preconditioner generation took: " << (g2o::get_time() - timing) * 1000 << "ms."<< endl;
  
      //CUDA stuff
      double *d_zm1, *d_zm2, *d_rm2;
      int *d_colIndsICP, *d_rowPtrsICP;
      double *d_valsICP; 
      double *cu_y;
  
      cutilSafeCall( cudaMalloc((void**)&cu_y, (_colum_count)*sizeof(double)) );  
      cutilSafeCall( cudaMalloc((void**)&d_colIndsICP, nzICP*sizeof(int)) );
      cutilSafeCall( cudaMalloc((void**)&d_valsICP, nzICP*sizeof(double)) );
      cutilSafeCall( cudaMalloc((void**)&d_rowPtrsICP, (_colum_count+1)*sizeof(int)) );
      cutilSafeCall( cudaMalloc((void**)&d_zm1, (_colum_count)*sizeof(double)) );
      cutilSafeCall( cudaMalloc((void**)&d_zm2, (_colum_count)*sizeof(double)) );
      cutilSafeCall( cudaMalloc((void**)&d_rm2, (_colum_count)*sizeof(double)) );
  
      cudaMemcpy(d_colIndsICP, colIndsICP, nzICP*sizeof(int), cudaMemcpyHostToDevice);
      cudaMemcpy(d_rowPtrsICP, rowPtrsICP, (_colum_count+1)*sizeof(int), cudaMemcpyHostToDevice);
      cudaMemcpy(d_valsICP, valsICP, nzICP*sizeof(double), cudaMemcpyHostToDevice);
  
        // create the analysis info object for the Non-Transpose case
        cusparseSolveAnalysisInfo_t info = 0;
        _status = cusparseCreateSolveAnalysisInfo(&info);
        if ( checkStatus ( _status, (char*)"cusparseCreateSolveAnalysisInfo" ) ) return EXIT_FAILURE;
    
        // create the analysis info object for the Transpose case
        cusparseSolveAnalysisInfo_t infoTrans = 0;
        _status = cusparseCreateSolveAnalysisInfo(&infoTrans);
        if ( checkStatus ( _status, (char*)"cusparseCreateSolveAnalysisInfo Trans" ) ) return EXIT_FAILURE;
    
        // Perform the analysis for the Non-Transpose case
        _status = cusparseDcsrsv_analysis(_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, _colum_count, _descr, d_valsICP, d_rowPtrsICP, d_colIndsICP, info);
        checkStatus ( _status, (char*)"susparseScsrv_analysis" );
    
        // Perform the analysis for the Transpose case
        _status = cusparseDcsrsv_analysis(_handle, CUSPARSE_OPERATION_TRANSPOSE, _colum_count, _descr, d_valsICP, d_rowPtrsICP, d_colIndsICP, infoTrans);
        if ( checkStatus ( _status, (char*)"susparseScsrv_analysis Trans" ) ) return EXIT_FAILURE;

      /* Preconditioned Conjugate Gradient.  
       Follows the description by Golub & Van Loan, "Matrix Computations 3rd ed.", Algorithm 10.3.1  */
      _iteration = 0;
      double r1 = cublasDdot(_colum_count, cu_b, 1, cu_b, 1);
      double r0 = 0;
      double alpha, beta;
      while (r1 > _tol*_tol && _iteration <= _max_iter) {
          // solve M z = H H^T z = r
  
          // Forward Solve: H y = r
          _status = cusparseDcsrsv_solve(_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, _colum_count, 1.0, _descr, 
                      d_valsICP, d_rowPtrsICP, d_colIndsICP, info, cu_b, cu_y);    
          if ( checkStatus ( _status, (char*)"susparseScsrv_solve" ) ) return EXIT_FAILURE;
          // Back Substitution: H^T z = y
          _status = cusparseDcsrsv_solve(_handle, CUSPARSE_OPERATION_TRANSPOSE, _colum_count, 1.0, _descr, 
                      d_valsICP, d_rowPtrsICP, d_colIndsICP, infoTrans, cu_y, d_zm1);    
          if ( checkStatus ( _status, (char*)"susparseScsrv_solve" ) ) return EXIT_FAILURE;
  
          _iteration++;
  
          if ( _iteration == 1 ) {
              cublasDcopy (_colum_count, d_zm1, 1, cu_p, 1);
          }
          else {
              beta = cublasDdot(_colum_count, cu_b, 1, d_zm1, 1)/cublasDdot(_colum_count, d_rm2, 1, d_zm2, 1);
              cublasDscal(_colum_count, beta, cu_p, 1);    
              cublasDaxpy (_colum_count, 1.0, d_zm1, 1, cu_p, 1) ;
          }
          _status = cusparseDcsrmv(_handle,CUSPARSE_OPERATION_NON_TRANSPOSE, _colum_count, _colum_count, 1.0, 
                          _descr, cu_A_ccs, cu_col_ptr_ccs, cu_row_idx_ccs, cu_p, 0.0, cu_Ap);
          if ( checkStatus ( _status, (char*)"cusparseDcsrmv" ) ) return EXIT_FAILURE;
          //cusparseDcsrmv(handle,CUSPARSE_OPERATION_NON_TRANSPOSE, _colum_count, _colum_count, 1.0, descr, d_val, cu_bow, d_col, cu_p, 0.0, cu_Ap);
          alpha = cublasDdot(_colum_count, cu_b, 1, d_zm1, 1)/cublasDdot (_colum_count, cu_p, 1, cu_Ap, 1);
          cublasDaxpy (_colum_count, alpha, cu_p, 1, cu_x, 1);
          cublasDcopy (_colum_count, cu_b, 1, d_rm2, 1);
          cublasDcopy (_colum_count, d_zm1, 1, d_zm2, 1);
          cublasDaxpy (_colum_count, -alpha, cu_Ap, 1, cu_b, 1);
          r1 = cublasDdot(_colum_count, cu_b, 1, cu_b, 1);
          shrLog("  iteration = %3d, residual = %e \n", _iteration, sqrt(r1));
      }

       shrLog("  iteration = %3d, residual = %e \n", _iteration, sqrt(r1));
       cerr << "gpu part took " << (get_time() - timing)*1000 << "ms." << endl;

       timing = get_time();
       cudaMemcpy(x, cu_x, sizeof(double)*_colum_count, cudaMemcpyDeviceToHost);
       cerr << "solution gpu -> main memory took " << (get_time() - timing)*1000 << "ms." << endl;

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





    protected:

  };


}// end namespace

#endif
