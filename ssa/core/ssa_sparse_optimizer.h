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

#ifndef __SSA_SPARSE_OPTIMIZER__
#define __SSA_SPARSE_OPTIMIZER__

#include "g2o/core/sparse_optimizer.h"
#include <QObject>

namespace ssa {

  struct SSASparseOptimizer :  public QObject, g2o::SparseOptimizer 
  {
    Q_OBJECT
    public:
    SSASparseOptimizer();
    ~SSASparseOptimizer();

    virtual void postIteration(int iteration); 

    signals:
       void iterationDone();
  };

 
} //end namespace

#endif
