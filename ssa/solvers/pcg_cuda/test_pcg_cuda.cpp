#include "linear_solver_pcg_cuda.h"
#include "ssa/core/ssa_graph_3d.h"
#include "ssa/core/allocate_solver.h"
#include "ssa/core/sparse_surface_adjustment.h"
#include "g2o/core/sparse_block_matrix.h"

using namespace std;
using namespace g2o;
using namespace Eigen;

typedef SparseBlockMatrix< MatrixXd > SparseBlockMatrixX;


std::ostream& operator << (std::ostream& os, const SparseBlockMatrixX::SparseMatrixBlock& m) {
  for (int i=0; i<m.rows(); ++i){
    for (int j=0; j<m.cols(); ++j)
      cerr << m(i,j) << " ";
    cerr << endl;
  }
  return os;
}

int main(int argc, const char* argv[])
{
 
  g2o::BlockSolverX::LinearSolverType* linearSolver = new g2o::LinearSolverCUDACG<g2o::BlockSolverX::PoseMatrixType>;
  //ssa.setSolver(linearSolver);

  //Eigen::Matrix4d matrix;
  //A = [ 1 7 0 0 ]
  //    [ 0 2 8 0 ]
  //    [ 5 0 3 9 ]
  //    [ 0 6 0 4 ]

//   int rcol[] = {1,2,3,4};
//   int ccol[] = {1,2,3,4};
//   cerr << "creation" << endl; 
//   SparseBlockMatrixX* M=new SparseBlockMatrixX(rcol, ccol, 4,4);

//   cerr << "block access" << endl;
// 
//   SparseBlockMatrixX::SparseMatrixBlock* b=M->block(0,0, true);
//   cerr << b->rows() << " " << b->cols() << endl;
//   for (int i=0; i<b->rows(); ++i){
//     for (int j=0; j<b->cols(); ++j){
//      (*b)(i,j) = 1;
//     }
//   }
// 
//   b=M->block(1,1, true);
//   (*b)(0,0) = 2;
// 
//   b=M->block(1,2, true);
//   (*b)(0,0) = 7;
// 
//   b=M->block(2,2, true);
//   (*b)(0,0) = 3;
// 
//   b=M->block(3,3, true);
//   (*b)(0,0) = 4;
// 
//   b=M->block(0,2, true);
//   (*b)(0,0) = 5;
// 
//   b=M->block(1,3, true);
//   (*b)(0,0) = 6;
  int rcol[] = {4};
  int ccol[] = {4};
  cerr << "creation" << endl;
  SparseBlockMatrixX* M=new SparseBlockMatrixX(rcol, ccol, 1,1);

 cerr << "block access" << endl;

  SparseBlockMatrixX::SparseMatrixBlock* b=M->block(0,0, true);
  cerr << b->rows() << " " << b->cols() << endl;
  (*b)(0,0)=3.0;
  (*b)(1,1)=4.0;
  (*b)(2,2)=10.0;
  (*b)(3,3)=3.0;

  (*b)(1,0)=1.0;
  (*b)(2,1)=1.0;
  (*b)(3,1)=3.0;

  (*b)(0,1)=1.0;
  (*b)(1,2)=1.0;
  (*b)(1,3)=3.0;


//   cerr << "block access 2" << endl;
//   b=M->block(0,2, true);
//   cerr << b->rows() << " " << b->cols() << endl;
//   for (int i=0; i<b->rows(); ++i)
//     for (int j=0; j<b->cols(); ++j){
//       (*b)(i,j)=i*b->cols()+j;
//     }
// 
//   b=M->block(3,2, true);
//   cerr << b->rows() << " " << b->cols() << endl;
//   for (int i=0; i<b->rows(); ++i)
//     for (int j=0; j<b->cols(); ++j){
//       (*b)(i,j)=i*b->cols()+j;
//     }

  std::cerr << (*M) << endl;


        int n = M->cols();
        int m = M->cols();

        MatrixXd H(n,m);
        H.setZero();


        int c_idx = 0;


        for (size_t i = 0; i < M->blockCols().size(); ++i)
        {
          int c_size = M->colsOfBlock(i);
          int r_idx = 0;

          const typename SparseBlockMatrix<MatrixXd>::IntBlockMap& col
              = M->blockCols()[i];
          if (col.size() > 0)
          {
            typename SparseBlockMatrix<MatrixXd>::IntBlockMap::const_iterator it;
            for (it = col.begin(); it != col.end(); ++it)
            {

              if (it->first <= (int)i)  // only the upper triangular block is needed
             {
                int r_size = M->rowsOfBlock(it->first);
                H.block(r_idx,c_idx,r_size,c_size)
                    = *(it->second);

                r_idx += r_size;
              }
            }
          }


          c_idx += c_size;
        }

  cerr << H << endl;

  int  _colum_count = M->cols();
  int  _non_zero_entries = M->nonZeros();
  int*  _row_indices_ptr= (int*)malloc(sizeof(int)*(_colum_count+1));                               // csr row pointers for matrix A
  int*  _colum_indices = (int*)malloc(sizeof(int)*_non_zero_entries);

  std::cerr << "_non_zero_entries: " << _non_zero_entries << std::endl;
  std::cerr << "_colum_count: " << _colum_count << std::endl;

   // copy Data from A to A_float
   double *A_ccs = NULL;
   A_ccs = (double*)malloc(sizeof(double)*_non_zero_entries);                            // csr values for matrix A
   _non_zero_entries = M->fillCCS(_row_indices_ptr, _colum_indices, A_ccs, false);
   _row_indices_ptr[_colum_count] = _non_zero_entries;

  cerr << "data =   \t";
  for(int i=0; i < _non_zero_entries; ++i)
    cerr << A_ccs[i] << " ";
  cerr << endl << "indices = \t";
  for(int i=0; i < _non_zero_entries; ++i)
    cerr << _colum_indices[i] << " ";
  cerr << endl << "pointer = \t";
  for(int i=0; i <= _colum_count; ++i)
    cerr << _row_indices_ptr[i] << " ";
  cerr << endl;

  double *x,*Mb;
  Mb = (double*)malloc(sizeof(double)*_colum_count); 
  x = (double*)malloc(sizeof(double)*_colum_count); 
  for(int i=0; i < _colum_count; ++i)
   Mb[i] = 1;

  linearSolver->solve(*M, x, Mb);
  for(int i=0; i < _colum_count; ++i)
    cerr << x[i] << ",";
    cerr << endl;
}
