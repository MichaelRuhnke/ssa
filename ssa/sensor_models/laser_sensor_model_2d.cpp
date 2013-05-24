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

#include "laser_sensor_model_2d.h"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace ssa{

  LaserSensorModel2D::LaserSensorModel2D(){};
  LaserSensorModel2D::~LaserSensorModel2D(){};

  void LaserSensorModel2D::applySensorModel(SparseSurfaceAdjustmentGraph2D& ssaGraph, LaserSensorParams& params)
  {
    for(unsigned int j = 0; j < ssaGraph._edges_observations.size(); ++j){
      EdgeSE2PointXYCov* edge = ssaGraph._edges_observations[j];
      double beamAngle = atan2(edge->measurement()(1), edge->measurement()(0));
      VertexPointXYCov* vertex = dynamic_cast<VertexPointXYCov*>(edge->vertices()[1]);
      setInformationForVertexPoint(edge, vertex, params);
    }
  }


  Eigen::Matrix2d
  LaserSensorModel2D::getCovarianceForPoint(VertexPointXYCov*& point, double& beamAngle, LaserSensorParams& params)
  {
    //calculate incidence angle in point frame
    Eigen::Vector2d beam = point->parentVertex()->estimate().translation() - point->estimate(); //TODO: Use sensor pose instead of robot pose
    Eigen::Vector2d normal = point->normal();
    beam = beam.normalized();
    normal = normal.normalized();
    double incidenceAngle = 0.0;
    if(point->covariance() != Eigen::Matrix2d::Identity())
           incidenceAngle = (acos(((normal.dot(beam)))));

    //calculate beam in parent vertex frame
    beam = point->parentVertex()->estimate().inverse() * point->estimate();
    return getCovarianceForPoint(beam, beamAngle, incidenceAngle, params);
  }

  void
  LaserSensorModel2D::setInformationForVertexPoint(EdgeSE2PointXYCov*& edge, VertexPointXYCov*& point, LaserSensorParams& params)
  {
    double beamAngle = atan2(edge->measurement()(1), edge->measurement()(0));
    //calculate incidence angle in point frame
    Eigen::Vector2d beam = edge->measurement();

    Eigen::Vector2d normal = point->normal();
    beam = beam.normalized();
    normal = normal.normalized();
    double incidenceAngle = 0.0;
    if(point->covariance() != Eigen::Matrix2d::Identity())
      incidenceAngle = acos(normal.dot(beam));

    double d = (tan(params.angularResolution * 0.5) * beam.norm());
    ///uncertainty of the surface measurement in direction of the beam
    double dk = d / cos(incidenceAngle);
    edge->information().setIdentity();
    edge->information()(0,0) = 1.0 / ( (dk + params.sensorPrecision) * (dk + params.sensorPrecision));
    double cError = 0.001 * beam.norm();
    edge->information()(1,1) = 1.0 / (cError * cError);

    Eigen::Rotation2Dd rot(beamAngle);
    Eigen::Matrix2d mrot = rot.toRotationMatrix();
    edge->information() = mrot * edge->information() * mrot.transpose();
  }

  Eigen::Matrix2d
  LaserSensorModel2D::getCovarianceForPoint(Eigen::Vector2d& beam, double& beamAngle, double& angle, LaserSensorParams& params)
  {
    double k11 = params.scale;
    //double k22 = params.scale;

    //diameter of the beam traversal region
    double d = (tan(params.angularResolution * 0.5) * beam.norm());
    //uncertainty of the surface measurement in direction of the beam
    double dk = d / cos(angle);
    //params.scale *
    double sigma11 = (k11 * (dk) * params.sensorPrecision);
    sigma11 *= sigma11;
    //params.scale * ;
    double cError = 0.0001;
    double sigma22 = (beam.norm() * cError);
    sigma22 *= sigma22;

    Eigen::Matrix2d covariance;
    covariance.fill(.0);
    covariance(0,0) = std::min(1.0, sigma11);
    covariance(1,1) = std::min(1.0, sigma22);

    //rotate into beam coordinate frame
    Eigen::Rotation2Dd rot(beamAngle);
    Eigen::Matrix2d mrot = rot.toRotationMatrix();
    covariance = mrot * covariance * mrot.transpose();
    return covariance;
  }


}