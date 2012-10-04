// g2o - General Graph Optimization
// Copyright (C) 2011 M. Ruhnke
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

#ifndef SOLVER_CUSP_CONJUGATE_GRADIENT_H
#define SOLVER_CUSP_CONJUGATE_GRADIENT_H

namespace g2o {
  using namespace std;
  /**
   * \brief linear solver using PCG, pre-conditioner is block Jacobi
   */
  class SolverCUSPCG 
  {
    public:
      SolverCUSPCG();
      ~SolverCUSPCG();

      //(A_ccs, _colum_ptr, _row_indices, _colum_count, _non_zero_entries, x , b)
      bool solve(double* A, int* ptr, int* indices, int cols, int nz, double* x, double* b);

      bool solveWithTiming(double* A, int* ptr, int* indices, int cols, int nz, double* x, double* b);
    protected:

  };


}// end namespace

#endif
