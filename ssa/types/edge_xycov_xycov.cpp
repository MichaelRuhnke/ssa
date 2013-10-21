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

#include "edge_xycov_xycov.h"

namespace ssa {
 using namespace Eigen;
 using namespace g2o;
 EdgePointXYCovPointXYCov::EdgePointXYCovPointXYCov() :
    BaseBinaryEdge<2, Eigen::Vector2d, VertexPointXYCov, VertexPointXYCov>()
  {
  }

  void EdgePointXYCovPointXYCov::computeError()
  {
    VertexPointXYCov* l1 = static_cast<VertexPointXYCov*>(_vertices[0]);
    VertexPointXYCov* l2 = static_cast<VertexPointXYCov*>(_vertices[1]);
    _error = l1->estimate() - l2->estimate();

    Vector2d x(1.0, 0.0);
    Vector2d normal = l1->globalNormal() + l2->globalNormal();
//     double diff = rad2deg(acos(l1->globalNormal().normalized().dot(l2->globalNormal().normalized())));
//     diff = 1.0 - (std::min(diff, 20.0 )/ 20.0);
    double norm = normal.norm();
    if(norm > 0){
      normal /= norm;
    } else {
      std::cerr << __PRETTY_FUNCTION__ << " Detectected undefined normals! " << PVAR(chi2()) << " " <<  normal << " " <<   " " << std::endl;
      information().setIdentity();
      return;
    }

    double angle = acos(normal.dot(x));
    information().setIdentity();
    information()(0,0) = 1000; /// corresponds to uncertainty of 1mm
    information()(1,1) = 1;
    Eigen::Rotation2Dd rot(-angle);
    information() = rot.matrix() * information() * rot.matrix().transpose();
    //information() = l1->covariance().inverse() + l2->covariance().inverse();
    if(chi2() < 0.0 || std::isnan(chi2())){
      std::cerr << PVAR(chi2()) << " " <<  information() << " " <<   " " << std::endl;
    }
  }

  bool EdgePointXYCovPointXYCov::read(std::istream& is)
  {
    Eigen::Vector2d p;
    is >> p[0] >> p[1];
    setMeasurement(p);

    is >> information()(0,0) >> information()(0,1) >> information()(1,1);
    information()(1,0) = information()(0,1);
    return true;
  }

  bool EdgePointXYCovPointXYCov::write(std::ostream& os) const
  {
    os << measurement()[0] << " " << measurement()[1] << " ";
    os << information()(0,0) << " " << information()(0,1) << " " << information()(1,1);
    return os.good();
  }

  void EdgePointXYCovPointXYCov::initialEstimate(const OptimizableGraph::VertexSet& from, OptimizableGraph::Vertex* to)
  {
    assert(from.size() == 1 && from.count(_vertices[0]) == 1 && "Can not initialize VertexSE2 position by VertexPointXYCov");
    VertexPointXYCov* vi = static_cast<VertexPointXYCov*>(_vertices[0]);
    VertexPointXYCov* vj = static_cast<VertexPointXYCov*>(_vertices[1]);
    if (from.count(vi) > 0 && to == vj) {
      vj->setEstimate(vi->estimate());
    }
  }

#ifndef NUMERIC_JACOBIAN_TWO_D_TYPES
  void EdgePointXYCovPointXYCov::linearizeOplus()
  {
//     const VertexPointXYCov* vi = static_cast<const VertexPointXYCov*>(_vertices[0]);
//     const VertexPointXYCov* vj = static_cast<const VertexPointXYCov*>(_vertices[1]);
//     const double& x1           = vi->estimate()[0];
//     const double& y1           = vi->estimate()[1];
//     const double& th1          = 0.0;
//     //vi->estimate().rotation().angle();
//     const double& x2           = vj->estimate()[0];
//     const double& y2           = vj->estimate()[1];
//
//     double aux_1 = cos(th1) ;
//     double aux_2 = -aux_1 ;
//     double aux_3 = sin(th1) ;

    _jacobianOplusXi( 0 , 0 ) = 1;
    _jacobianOplusXi( 0 , 1 ) = 0;
    _jacobianOplusXi( 1 , 0 ) = 0;
    _jacobianOplusXi( 1 , 1 ) = 1;

    _jacobianOplusXj( 0 , 0 ) = -1;
    _jacobianOplusXj( 0 , 1 ) = 0;
    _jacobianOplusXj( 1 , 0 ) = 0;
    _jacobianOplusXj( 1 , 1 ) = -1;
  }
#endif

} //end namespace