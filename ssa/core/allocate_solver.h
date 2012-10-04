// Sparse Surface Optimization
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
// 
// SSA is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SSA is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __SSA_ALLOCATE_SOLVER__
#define __SSA_ALLOCATE_SOLVER__

#include "g2o/solvers/pcg/linear_solver_pcg.h"
#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/csparse/linear_solver_csparse.h"
//#include "g2o/solvers/cusp/linear_solver_cuda_conjugate_gradient.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

  template<typename T>
  typename T::LinearSolverType* AllocateLinearSolver(int solver = 0)
  {
    switch(solver){
      case 1:
        //CHOLMOD
        std::cerr << "Allocate CHOLMOD solver..." << std::endl;
        return new g2o::LinearSolverCholmod<typename T::PoseMatrixType>;
        break;
      case 2:
        //PCG
        std::cerr << "Allocate PCG solver..." << std::endl;
        return new g2o::LinearSolverPCG<typename T::PoseMatrixType>;
        break;
//       case 3:
//         //CUDA
//         std::cerr << "Allocate CUDA Conjugate Gradient solver..." << std::endl;
//         return new g2o::LinearSolverCUDACG<typename T::PoseMatrixType>;
//         break;
      default:
        std::cerr << "Allocate CSPARSE solver..." << std::endl;
        return new g2o::LinearSolverCSparse<typename T::PoseMatrixType>;
    }
  }

}
#endif
