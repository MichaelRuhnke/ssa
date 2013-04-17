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

#include "edge_xyzcov_xyzcov.h"

namespace ssa {
 using namespace Eigen;
 using namespace g2o;
 EdgePointXYZCovPointXYZCov::EdgePointXYZCovPointXYZCov() :
    BaseBinaryEdge<3, Eigen::Vector3d, VertexPointXYZCov, VertexPointXYZCov>()
  {
  }

  void EdgePointXYZCovPointXYZCov::computeError()
  {
    const VertexPointXYZCov* l1 = static_cast<const VertexPointXYZCov*>(_vertices[0]);
    const VertexPointXYZCov* l2 = static_cast<const VertexPointXYZCov*>(_vertices[1]);
    _error = l1->estimate() - l2->estimate();

    /** this is the way Kurt computes the information matrix
      // re-define the information matrix
      // topLeftCorner<3,3>() is the rotation()
      const Matrix3d transform = ( vp0->estimate().inverse() *  vp1->estimate() ).matrix().topLeftCorner<3,3>();
      information() = ( cov0 + transform * cov1 * transform.transpose() ).inverse();
    **/  
//     Matrix3d transform = ( l1->parentVertex()->estimate().inverse() *  l2->parentVertex()->estimate() ).matrix().topLeftCorner<3,3>();
//     information() = ( l1->covariance() + l2->covariance() ).inverse(); ///THIS SHOULD BE RIGHT after changing covariances to local frame
    information() = (l1->covariance().inverse() + l2->covariance().inverse()) ;

    if(chi2() < 0.0){
      std::cerr << PVAR(chi2()) << " " <<  l1->covariance().inverse() << " " <<  l2->covariance().inverse() << " " << std::endl;
    }
  }

  bool EdgePointXYZCovPointXYZCov::read(std::istream& is)
  {
    /** the measurement should always be 0,0,0 and the information matrix is computed in computeError() */
    Eigen::Vector3d p(0.0,0.0,0.0);
//     is >> p[0] >> p[1] >> p[2];
    setMeasurement(p);
// 
//     for (int i=0; i<3; i++)
//       for (int j=i; j<3; j++) {
//         is >> information()(i,j);
//         if (i!=j)
//           information()(j,i)=information()(i,j);
//       }
    return true;
  }

  bool EdgePointXYZCovPointXYZCov::write(std::ostream& os) const
  {
    /** the measurement should always be 0,0,0 and the information matrix is computed in computeError() */
//     for (int i=0; i<3; i++)
//       os << measurement()[i] << " ";
//     for (int i=0; i<3; i++)
//       for (int j=i; j<3; j++){
//         os << " " <<  information()(i,j);
//       }
    return os.good();
  }

  void EdgePointXYZCovPointXYZCov::initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to)
  {
    assert(from.size() == 1 && from.count(_vertices[0]) == 1 && "Can not initialize VertexSE3 position by VertexPointXYZCov");
    VertexPointXYZCov* vi = static_cast<VertexPointXYZCov*>(_vertices[0]);
    VertexPointXYZCov* vj = static_cast<VertexPointXYZCov*>(_vertices[1]);
    if (from.count(vi) > 0 && to == vj) {
      vj->setEstimate(vi->estimate());
    }
  }

#ifndef NUMERIC_JACOBIAN_THREE_D_TYPES
  void EdgePointXYZCovPointXYZCov::linearizeOplus()
  {
    _jacobianOplusXi( 0 , 0 ) = 1;
    _jacobianOplusXi( 0 , 1 ) = 0;
    _jacobianOplusXi( 0 , 2 ) = 0;
    _jacobianOplusXi( 1 , 0 ) = 0;
    _jacobianOplusXi( 1 , 1 ) = 1;
    _jacobianOplusXi( 1 , 2 ) = 0;
    _jacobianOplusXi( 2,  0 ) = 0;
    _jacobianOplusXi( 2 , 1 ) = 0;
    _jacobianOplusXi( 2 , 2 ) = 1;

    _jacobianOplusXj( 0 , 0 ) = -1;
    _jacobianOplusXj( 0 , 1 ) = 0;
    _jacobianOplusXj( 0 , 2 ) = 0;
    _jacobianOplusXj( 1 , 0 ) = 0;
    _jacobianOplusXj( 1 , 1 ) = -1;
    _jacobianOplusXj( 1 , 2 ) = 0;
    _jacobianOplusXj( 2 , 0 ) = 0;
    _jacobianOplusXj( 2 , 1 ) = 0;
    _jacobianOplusXj( 2 , 2 ) = -1;
  }
#endif

} //end namespace