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

//#include <Eigen/Eigenvalues>
#include "g2o/types/slam3d/types_slam3d.h"
#include "pcl/common/eigen.h"
#include "vertex_point_xyzcov.h"
#include "edge_se3_xyzcov.h"

namespace ssa {
  using namespace Eigen;
  using namespace g2o;

  VertexPointXYZCov::VertexPointXYZCov() : VertexPointXYZ(), 
    _parentVertex(0), 
    _cov(Eigen::Matrix3d::Identity()), 
    _normal(Eigen::Vector3d::Zero()), 
    _hasNormal(false)
  {
    
  };

  VertexPointXYZCov::VertexPointXYZCov(int id, Eigen::Vector3d estimate) : VertexPointXYZ(),  
    _parentVertex(0),
     _cov(Eigen::Matrix3d::Identity()), 
    _normal(Eigen::Vector3d::Zero()), 
    _hasNormal(false)
  {
    _id = id;
    _estimate = estimate;
  }

  bool VertexPointXYZCov::read(std::istream& is)
  {
    ///read id of sensor pose vertex
    is >> _parentVertexId;
    
    ///read point pose
    Eigen::Vector3d p;
    is >> p[0] >> p[1] >> p[2];
    setEstimate(p);
    
    ///read color information
    uint r,g,b;
    is >> r >> g >> b; 
    cr = (unsigned char) r;
    cg = (unsigned char) g;
    cb = (unsigned char) b;
    
    ///read point covariance
    for (int i=0; i<3; i++)
      for (int j=i; j<3; j++) {
        is >> _cov(i,j);
        if (i!=j)
          _cov(j,i)=_cov(i,j);
      }
    return true;
  }

  bool VertexPointXYZCov::write(std::ostream& os) const
  {
    /// we store the id of the observation / parent vertex for faster lookups... 
    /// this information is also available through the connected edges 
    os << _parentVertex->id() << " ";
    
    /// estimated pose of point (global coordinate frame): 
    /// this will be different from the measured depth after optimization!
    for (int i=0; i<3; i++)
      os << estimate()[i] << " ";

    /// color information for the point (if available) 
    uint r,g,b;
    r = (uint) cr;
    g = (uint) cg;
    b = (uint) cb;
    os << r << " " << g << " " << b << " ";   
    
    ///covariance of the point, computed from the neighborhood (do we really need to store in on the drive??)
    for (int i=0; i<3; i++)
      for (int j=i; j<3; j++){
        os << _cov(i,j) << " ";
      }
      
    return os.good();
  }

    g2o::VertexSE3* VertexPointXYZCov::parentVertex() {
    if(_parentVertex != 0)
      return _parentVertex;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=edges().begin(); it!=edges().end(); it++){
      EdgeSE3PointXYZCov* e2=dynamic_cast<EdgeSE3PointXYZCov*>(*it);
      if(e2){
        g2o::VertexSE3* p = dynamic_cast<g2o::VertexSE3* >(e2->vertices()[0]);
	if(p){
	  _parentVertex = p;
	}
        return p;
      }
    }
    g2o::VertexSE3* p = 0;
    return p;
  }
  
  g2o::VertexSE3* VertexPointXYZCov::parentVertex() const{
    if(_parentVertex != 0)
      return _parentVertex;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=edges().begin(); it!=edges().end(); it++){
      EdgeSE3PointXYZCov* e2=dynamic_cast<EdgeSE3PointXYZCov*>(*it);
      if(e2){
        g2o::VertexSE3* p = dynamic_cast<g2o::VertexSE3* >(e2->vertices()[0]);
        return p;
      }
    }
    g2o::VertexSE3* p = 0;
    return p;
  }

  std::vector< g2o::VertexSE3* > VertexPointXYZCov::parentVertices() const{
    std::vector< g2o::VertexSE3* > parents;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=edges().begin(); it!=edges().end(); it++){
      EdgeSE3PointXYZCov* e2=dynamic_cast<EdgeSE3PointXYZCov*>(*it);
      if(e2){
        g2o::VertexSE3* p = dynamic_cast<g2o::VertexSE3* >(e2->vertices()[0]);
        if(p)
          parents.push_back(p);
      }
    }
    return parents;
  }

  void VertexPointXYZCov::setParentVertex(VertexSE3* pose){
    _parentVertex = pose;
  }

  void VertexPointXYZCov::updateNormal(Eigen::Matrix3d& cov){
    Eigen::Matrix3d eigvectors;
    Eigen::Vector3d eigvalues;
    if(cov != Eigen::Matrix3d::Identity()){
      pcl::eigen33<Matrix3d, Vector3d> (cov, eigvectors, eigvalues);

      _normal = Vector3d(0.0, 0.0, 0.0);
      double eig_sum = eigvalues.sum();
      if(eig_sum != 0){
        if(eigvalues(0) < eigvalues(1) && eigvalues(0) < eigvalues(2)){
          _normal = Vector3d(eigvectors(0,0), eigvectors(1,0), eigvectors(2,0));
        }
        if(eigvalues(1) < eigvalues(0) && eigvalues(1) < eigvalues(2)){
          _normal = Vector3d(eigvectors(0,1), eigvectors(1,1), eigvectors(2,1));
        }
        if(eigvalues(2) < eigvalues(0) && eigvalues(2) < eigvalues(1)){
          _normal = Vector3d(eigvectors(0,2), eigvectors(1,2), eigvectors(2,2));
        }
  
        //check direction of normal
        Eigen::Vector3d laserPoint(_estimate[0], _estimate[1], _estimate[2]);
        Eigen::Vector3d beam = (parentVertex()->estimate().translation() - laserPoint);
        if(beam.dot(_normal) < 0){
          _normal = -_normal;
        }
        
//         Eigen::Matrix3d coordinateFrameRotation;
//         Eigen::Vector3d yDirection(0,1,0);
// 	coordinateFrameRotation.row(2) = _normal;
// 	yDirection = yDirection - _normal(1)*_normal;
// 	yDirection.normalize();/// need to check if y is close to 0
// 	coordinateFrameRotation.row(1) = yDirection;
// 	coordinateFrameRotation.row(0) = _normal.cross(yDirection);
// 	
// 	_cov = Matrix3d::Identity();
// 	_cov(0,0) = 0.0001; ///Maybe we should make this resolution dependent
// 	_cov(1,1) = 0.0001; ///Maybe we should make this resolution dependent
// 	_cov(2,2) = 0.0001; ///Maybe we should make this resolution dependent
// 	_cov = coordinateFrameRotation.transpose()*_cov*coordinateFrameRotation;	
	_normal = parentVertex()->estimate().linear().inverse() * _normal; //rotate in local coordinate frame
      }     
      _hasNormal = true;
    }  
  }
  
  Eigen::Vector3d& VertexPointXYZCov::normal(){
    Eigen::Matrix3d eigvectors;
    Eigen::Vector3d eigvalues;
    if(!_hasNormal && _cov != Eigen::Matrix3d::Identity()){
      pcl::eigen33<Matrix3d, Vector3d> (_cov, eigvectors, eigvalues);

      _normal = Vector3d(0.0, 0.0, 0.0);
      double eig_sum = eigvalues.sum();
      if(eig_sum != 0){
        if(eigvalues(0) < eigvalues(1) && eigvalues(0) < eigvalues(2)){
          _normal = Vector3d(eigvectors(0,0), eigvectors(1,0), eigvectors(2,0));
        }
        if(eigvalues(1) < eigvalues(0) && eigvalues(1) < eigvalues(2)){
          _normal = Vector3d(eigvectors(0,1), eigvectors(1,1), eigvectors(2,1));
        }
        if(eigvalues(2) < eigvalues(0) && eigvalues(2) < eigvalues(1)){
          _normal = Vector3d(eigvectors(0,2), eigvectors(1,2), eigvectors(2,2));
        }
  
        //check direction of normal
        Vector3d laserPoint(_estimate[0], _estimate[1], _estimate[2]);
        Vector3d beam = (parentVertex()->estimate().translation() - laserPoint);
        if(beam.dot(_normal) < 0){
          _normal = -_normal;
        }
      }     
      _hasNormal = true;
    }
    return _normal;
  }
  
  Eigen::Vector3d VertexPointXYZCov::globalNormal(){
    Vector3d normal(0.0, 0.0, 0.0);
    if(_hasNormal)
	normal = parentVertex()->estimate().linear() * _normal;
    return normal;
  }

  Vector3d VertexPointXYZCov::normal() const{
    return _normal;
  }

  Eigen::Matrix3d VertexPointXYZCov::covariance() const{
    return _cov;
  }

  Eigen::Matrix3d& VertexPointXYZCov::covariance(){
    return _cov;
  }

  Eigen::Matrix4d VertexPointXYZCov::covTransform(){
    EigenSolver<Matrix3d> solv(_cov);
    Matrix3d eigenVectors = solv.eigenvectors().real();

    Eigen::Matrix4d covGlTransform;
    covGlTransform.setIdentity();
    for(int i=0; i < 3; ++i)
      for(int j=0; j < 3; ++j)
        covGlTransform(i,j) = eigenVectors(i,j);

    for(int i=0; i < 3; ++i)
        covGlTransform(i,3) = estimate()(i);

    return covGlTransform; 
  }
  
} //end namespace

