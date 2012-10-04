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

#ifndef LINEAR_SOLVER_CUSP_CONJUGATE_GRADIENT_H
#define LINEAR_SOLVER_CUSP_CONJUGATE_GRADIENT_H

#include "g2o/core/linear_solver.h"
#include "g2o/core/batch_stats.h"

#include <vector>
#include "pcg_cusp.h"

#include "EXTERNAL/g2o/g2o/stuff/timeutil.h"

namespace g2o {
  using namespace std;
  /**
   * \brief linear solver using PCG, pre-conditioner is block Jacobi
   */
  template <typename MatrixType>
  class LinearSolverCuspCG : public LinearSolver<MatrixType>
  {
    public:
      LinearSolverCuspCG() :
      LinearSolver<MatrixType>()
      {
      }

      virtual ~LinearSolverCuspCG()
      {
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


      protected:
        int  _colum_count;
        int  _non_zero_entries;
        int* _colum_ptr;
        int* _row_indices;
        int _max_iter ;
        int _iteration;
        double _tol;

      public:



      bool solve(const SparseBlockMatrix<MatrixType>& A, double* x, double* b)
      {
        double timing = get_time();
        /**  initialize ccs colum pointers for matrix A */
        if(_colum_count < A.cols()){
          delete _colum_ptr;
          _colum_ptr = (int*)malloc(sizeof(int)*(A.cols()+1));
        }
        _colum_count = A.cols();

        int nonZeros = A.nonZeros() * 2; //assign double of the non zero space for upper to full ccs
        /** initialize ccs row indices for matrix A */
        if(_non_zero_entries < (int) nonZeros){
          delete _row_indices;
          _row_indices = (int*)malloc(sizeof(int)*(nonZeros));
        }
        _non_zero_entries = nonZeros;
        cerr << "memory allocation took " << (get_time() - timing) * 1000 << " ms." << endl;

        timing = get_time();
        // get data in CCS format
        double *A_ccs = NULL;
        A_ccs = (double*)malloc(sizeof(double)*nonZeros);
        A.fillCCS(_colum_ptr, _row_indices, A_ccs, true);
        ccsUpperToccsFull(A_ccs, _colum_ptr, _row_indices, _colum_count);

//         int n = A.cols();
//         int m = A.cols();
//         cerr << "cols" << A.cols() << " non_zeros" << _non_zero_entries << endl;
//         MatrixXd H(n,m);
//         H.setZero();
// 
//         int col = 0;
//         for(int i=0;i < A.cols(); ++i){
//           int c_start = _colum_ptr[i];
//           int c_end = _colum_ptr[i+1];
//           cerr << i << " " << c_start << " -> " << c_end << endl;
//           for(int j=c_start;j < c_end; ++j){
//             int row = _row_indices[j];
//             cerr << col << "," << row << endl;
//             H(col,row) = A_ccs[j];
//           }
//           col++;
//         }
// 
//         for(int i=0;i < n; ++i){
//           for(int j=0;j < n; ++j){
//             cout << H(i,j) << " ";
//           }
//           cout << endl;
//         }

        cerr << "fillCCS took " << (get_time() - timing) * 1000 << " ms." << endl;
//         exit(0);
        SolverCUSPCG cusp_cg;
        cusp_cg.solve(A_ccs, _colum_ptr, _row_indices, _colum_count, _non_zero_entries, x , b);

        delete A_ccs;
        return true;
      }

      int ccsUpperToccsFull(double* A_ccs, int* ptr_ccs, int* indices_ccs, int cols){
        double t = get_time();
        std::map<int,std::map<int, double> > a_coo;
        int col = 0;
        for(int i=0;i < cols; ++i){
          int c_start = ptr_ccs[i];
          int c_end = ptr_ccs[i+1];
          for(int j=c_start;j < c_end; ++j){
            int row = indices_ccs[j];
            a_coo[col][row] = A_ccs[j];
            if(col != row){
              a_coo[row][col] = A_ccs[j];
            }
          }
          col++;
        }
        cerr << "ccs -> std::map " << (get_time() - t) * 1000 << endl;
        int col_index = 0;
        int nz_index = 0;
        for(std::map<int,std::map<int, double> >::iterator it=a_coo.begin(); it!=a_coo.end(); ++it){
          ptr_ccs[col_index] = nz_index;
          for(std::map<int, double>::iterator itt=it->second.begin(); itt!=it->second.end(); ++itt){
            const int& j = itt->first;
            A_ccs[nz_index] = itt->second;
            indices_ccs[nz_index] = j;
            nz_index++;

          }
          col_index++;
        }
        ptr_ccs[col_index] = nz_index;
        cerr << "ccs upper to full ccs took " << (get_time() - t) * 1000 << " ms." << endl;
        return nz_index;
      }

    protected:

  };


}// end namespace

#endif
