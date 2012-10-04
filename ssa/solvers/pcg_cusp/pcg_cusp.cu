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

#include "pcg_cusp.h"

#include <cusp/precond/diagonal.h>
#include <cusp/transpose.h>
#include <cusp/elementwise.h>
#include <cusp/csr_matrix.h>
#include <cusp/krylov/cg.h>
#include <cusp/krylov/bicgstab.h>
 
#include "EXTERNAL/g2o/g2o/stuff/timeutil.h"

namespace g2o {

  SolverCUSPCG::SolverCUSPCG() { };
  SolverCUSPCG::~SolverCUSPCG() { };

  bool SolverCUSPCG::solve(double* A, int* ptr, int* indices, int cols, int nz, double* x, double* b){

    cusp::array1d<double,cusp::host_memory> x_cusp_host;
    cusp::csr_matrix<int,double,cusp::device_memory>  A_ccs_cusp_device(cols, cols, nz);
    cusp::array1d<double,cusp::device_memory> x_cusp_device(A_ccs_cusp_device.num_rows, 0.0);
    cusp::array1d<double,cusp::device_memory> b_cusp_device(A_ccs_cusp_device.num_rows);

    thrust::copy(A, A + nz, A_ccs_cusp_device.values.begin());
    thrust::copy(indices, indices + nz, A_ccs_cusp_device.column_indices.begin());
    thrust::copy(ptr, ptr + cols + 1, A_ccs_cusp_device.row_offsets.begin());
    thrust::copy(b, b + cols, b_cusp_device.begin());

//     cusp::csr_matrix<int,double,cusp::device_memory>  At_ccs_cusp_device(cols, cols, nz);
//     cusp::transpose(A_ccs_cusp_device, At_ccs_cusp_device);
//     cusp::csr_matrix<int,double,cusp::host_memory>  At_ccs_cusp_host;
//   
//     //removing diagonal of transposed matrix
//     At_ccs_cusp_host = At_ccs_cusp_device;
//     int c=0;
//     for(int i=0;i < (At_ccs_cusp_host.row_offsets.size()-1); ++i){
//       int c_start = At_ccs_cusp_host.row_offsets[i];
//       int c_end = At_ccs_cusp_host.row_offsets[i+1];
//       for(int j=c_start;j < c_end; ++j){
//         if(At_ccs_cusp_host.column_indices[j] == c)
//           At_ccs_cusp_host.values[j] = 0.0;
//       }
//       c++;
//     }
//     At_ccs_cusp_device = At_ccs_cusp_host;
//     cusp::add(A_ccs_cusp_device, At_ccs_cusp_device, A_ccs_cusp_device);

    cusp::convergence_monitor<double> monitor(b_cusp_device, cols, 1e-8);

    // set preconditioner (identity)
    //cusp::identity_operator<double,cusp::device_memory> M(A_ccs_cusp_device.num_rows, A_ccs_cusp_device.num_rows);

    // setup preconditioner Jacoby
    cusp::precond::diagonal<double, cusp::device_memory> M(A_ccs_cusp_device);

    // solve the linear system A * x = b with the Conjugate Gradient method
    cusp::krylov::cg(A_ccs_cusp_device, x_cusp_device, b_cusp_device, monitor, M);
    x_cusp_host = x_cusp_device;
    // copy x back
    for (int i = 0; i < cols; i++) {
      x[i] = x_cusp_host[i];
    } 
    return true;
  }

  bool SolverCUSPCG::solveWithTiming(double* A, int* ptr, int* indices, int cols, int nz, double* x, double* b){

    double timing = get_time();
    cerr << "allocating memory...  \t ";
    cusp::array1d<double,cusp::host_memory> x_cusp_host;

    // copy to the device
    cusp::csr_matrix<int,double,cusp::device_memory>  A_ccs_cusp_device(cols, cols, nz);
    cusp::array1d<double,cusp::device_memory> x_cusp_device(A_ccs_cusp_device.num_rows, 0.0);
    cusp::array1d<double,cusp::device_memory> b_cusp_device(A_ccs_cusp_device.num_rows);
    cerr << "done in " << (get_time() - timing) * 1000 << " ms." << endl;

    timing = get_time();
    cerr << "copy data to gpu...  \t";
    thrust::copy(A, A + nz, A_ccs_cusp_device.values.begin());
    thrust::copy(indices, indices + nz, A_ccs_cusp_device.column_indices.begin());
    thrust::copy(ptr, ptr + cols + 1, A_ccs_cusp_device.row_offsets.begin());
    thrust::copy(b, b + cols, b_cusp_device.begin());
    cerr << "done in " << (get_time() - timing) * 1000 << " ms." << endl;

//     timing = get_time();
//     cerr << "csr upper to csr full ... \t";
//     cusp::csr_matrix<int,double,cusp::device_memory>  At_ccs_cusp_device(cols, cols, nz);
//     //removing diagonal of transposed matrix
//     //At_ccs_cusp_host = At_ccs_cusp_device;
//     int col = 0;
//     for(int i=0;i < cols; ++i){
//       int c_start = ptr[i];
//       int c_end = ptr[i+1];
//       for(int j=c_start;j < c_end; ++j){
//         int row = indices[j];
//         if(col == row)
//           A[j] = 0.0;
//       }
//       col++;
//     }
//     thrust::copy(A, A + nz, At_ccs_cusp_device.values.begin());
//     At_ccs_cusp_device.column_indices = A_ccs_cusp_device.column_indices;
//     At_ccs_cusp_device.row_offsets = A_ccs_cusp_device.row_offsets;
// 
//     cusp::csr_matrix<int,double,cusp::device_memory>  Atmp(cols, cols, nz);
//     cusp::transpose(At_ccs_cusp_device, Atmp);
//     cusp::add(A_ccs_cusp_device, Atmp, A_ccs_cusp_device);
//     cerr << "done in " << (get_time() - timing) * 1000 << " ms." << endl;

    cerr << "running conjugate gradient with jacoby preconditioner... \t";
    // set stopping criteria:
    //  iteration_limit    = 100
    //  relative_tolerance = 1e-3
    cusp::convergence_monitor<double> monitor(b_cusp_device, cols, 1e-8);

    // set preconditioner (identity)
    //cusp::identity_operator<double,cusp::device_memory> M(A_ccs_cusp_device.num_rows, A_ccs_cusp_device.num_rows);

    // setup preconditioner Jacoby
    cusp::precond::diagonal<double, cusp::device_memory> M(A_ccs_cusp_device);

    // solve the linear system A * x = b with the Conjugate Gradient method
    cusp::krylov::cg(A_ccs_cusp_device, x_cusp_device, b_cusp_device, monitor, M);
    cerr << "done in " << (get_time() - timing) * 1000 << " ms." << endl;

    timing = get_time();
    cerr << "copy result back to host memory... \t";
    x_cusp_host = x_cusp_device;
    // copy x back
    for (int i = 0; i < cols; i++) {
      x[i] = x_cusp_host[i];
    } 
    cerr << "done in " << (get_time() - timing) * 1000 << " ms." << endl;
    return true;
  }

}// end namespace


