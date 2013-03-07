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

#include "vertex_point_xycov.h"
#include <Eigen/Eigenvalues>
#include "g2o/types/slam2d/vertex_se2.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

  VertexPointXYCov::VertexPointXYCov()
  {
    _parentVertex = 0;
    ratio = 0.0;
  };

  bool VertexPointXYCov::read(std::istream& is)
  {
    is >> _estimate[0] >> _estimate[1] >> _cov(0,0) >> _cov(0,1) >> _cov(1,1) >> _parentVertexId >> ratio;
    _cov(1,0) = _cov(0,1);
    return true;
  }

  bool VertexPointXYCov::write(std::ostream& os) const
  {
    os << _estimate[0] << " " << _estimate[1] << " " << _cov(0,0) << " " << _cov(0,1) << " " << _cov(1,1) << " " << _parentVertex->id() << " " << ratio;
    return os.good();
  }

  g2o::VertexSE2* VertexPointXYCov::parentVertex() const{
    return _parentVertex;
  }

  unsigned int VertexPointXYCov::parentVertexId() const{
    if(_parentVertex != 0)
      return _parentVertex->id();
    return _parentVertexId;
  }

  void VertexPointXYCov::setParentVertex(VertexSE2* pose){
    _parentVertex = pose;
  }

  void VertexPointXYCov::updateNormal(Eigen::Matrix2d& cov){
    Matrix2d eigvectors;
    Vector2d eigvalues;
    Eigen::EigenSolver<Matrix2d> solv(cov);
    eigvectors = solv.eigenvectors().real();
    eigvalues = solv.eigenvalues().real();
    
    if(eigvalues(1) <= eigvalues(0)){
      _normal = Vector2d(eigvectors(0,1), eigvectors(1,1));
    } else {
      _normal = Vector2d(eigvectors(0,0), eigvectors(1,0));
    }
    
    ///check direction of normal
    Vector2d laserPoint(_estimate[0], _estimate[1]);
    Vector2d robotPose(parentVertex()->estimate()[0], parentVertex()->estimate()[1]);

    if((laserPoint - robotPose).normalized().dot(_normal) > 0){
      _normal = -_normal;
    }
  }
      
  Vector2d& VertexPointXYCov::normal(){
    Matrix2d eigvectors;
    Vector2d eigvalues;
    if(!_hasNormal){

      //compute eigenvectors and eigenvalues
      Eigen::EigenSolver<Matrix2d> solv(_cov);
      eigvectors = solv.eigenvectors().real();
      eigvalues = solv.eigenvalues().real();

      if(eigvalues(1) <= eigvalues(0)){
        _normal = Vector2d(eigvectors(0,1), eigvectors(1,1));
      } else {
        _normal = Vector2d(eigvectors(0,0), eigvectors(1,0));
      }

      //check direction of normal
      Vector2d laserPoint(_estimate[0], _estimate[1]);
      Vector2d robotPose(parentVertex()->estimate()[0], parentVertex()->estimate()[1]);
      //Eigen::Rotation2Dd rot(parentVertex()->estimate()[2]); rot *
      if((laserPoint - robotPose).normalized().dot(_normal) > 0){
        _normal = -_normal;
      }
      Eigen::Rotation2Dd rot(parentVertex()->estimate()[2]);
      _normal = rot * _normal;
      _hasNormal = true;
    }
    return _normal;
  }

  Vector2d VertexPointXYCov::normal() const{
    return _normal;
  }

   Eigen::Vector2d VertexPointXYCov::globalNormal(){
    Vector2d globalNormal(0.0, 0.0);
    if(parentVertex()){
      Eigen::Rotation2Dd rot(parentVertex()->estimate()[2]);
      globalNormal = rot * normal();
    }
    return globalNormal;
  }

   Eigen::Matrix2d VertexPointXYCov::covariance() const{
    return _cov;
  }

   Eigen::Matrix2d& VertexPointXYCov::covariance(){
    return _cov;
  }
}

