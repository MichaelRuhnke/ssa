// Sparse Surface Optimization 2D
// Copyright (C) 2011 M. Ruhnke, R. Kuemmerle, G. Grisetti, W. Burgard
// 
// SSA2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published
// by the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// SSA2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "edge_se2_xycov.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

  EdgeSE2PointXYCov::EdgeSE2PointXYCov() :
    BaseBinaryEdge<2, Vector2d, VertexSE2, VertexPointXYCov>()
  {
  }

  bool EdgeSE2PointXYCov::read(std::istream& is)
  {
    Eigen::Vector2d p;
    is >> p[0] >> p[1];
    setMeasurement(p);

    is >> information()(0,0) >> information()(0,1) >> information()(1,1);
    information()(1,0) = information()(0,1);
    return true;
  }

  bool EdgeSE2PointXYCov::write(std::ostream& os) const
  {
    os << measurement()[0] << " " << measurement()[1] << " ";
    os << information()(0,0) << " " << information()(0,1) << " " << information()(1,1);
    return os.good();
  }

  void EdgeSE2PointXYCov::initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to)
  {
    assert(from.size() == 1 && from.count(_vertices[0]) == 1 && "Can not initialize VertexSE2 position by VertexPointXYCovCov");

    VertexSE2* vi     = static_cast<VertexSE2*>(_vertices[0]);
    VertexPointXYCov* vj = static_cast<VertexPointXYCov*>(_vertices[1]);
    if (from.count(vi) > 0 && to == vj) {
      vj->setEstimate(vi->estimate() * measurement());
    }
  }

#ifndef NUMERIC_JACOBIAN_TWO_D_TYPES
  void EdgeSE2PointXYCov::linearizeOplus()
  {
    const VertexSE2* vi     = static_cast<const VertexSE2*>(_vertices[0]);
    const VertexPointXYCov* vj = static_cast<const VertexPointXYCov*>(_vertices[1]);
    const double& x1        = vi->estimate().translation()[0];
    const double& y1        = vi->estimate().translation()[1];
    const double& th1       = vi->estimate().rotation().angle();
    const double& x2        = vj->estimate()[0];
    const double& y2        = vj->estimate()[1];

    double aux_1 = cos(th1) ;
    double aux_2 = -aux_1 ;
    double aux_3 = sin(th1) ;

    _jacobianOplusXi( 0 , 0 ) = aux_2 ;
    _jacobianOplusXi( 0 , 1 ) = -aux_3 ;
    _jacobianOplusXi( 0 , 2 ) = aux_1*y2-aux_1*y1-aux_3*x2+aux_3*x1 ;
    _jacobianOplusXi( 1 , 0 ) = aux_3 ;
    _jacobianOplusXi( 1 , 1 ) = aux_2 ;
    _jacobianOplusXi( 1 , 2 ) = -aux_3*y2+aux_3*y1-aux_1*x2+aux_1*x1 ;

    _jacobianOplusXj( 0 , 0 ) = aux_1 ;
    _jacobianOplusXj( 0 , 1 ) = aux_3 ;
    _jacobianOplusXj( 1 , 0 ) = -aux_3 ;
    _jacobianOplusXj( 1 , 1 ) = aux_1 ;
  }
#endif

}

