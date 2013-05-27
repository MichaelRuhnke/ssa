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
//       double beamAngle = atan2(edge->measurement()(1), edge->measurement()(0));
      VertexPointXYCov* vertex = dynamic_cast<VertexPointXYCov*>(edge->vertices()[1]);
      setInformationForVertexPoint(edge, vertex, params);
    }
  }

  void
  LaserSensorModel2D::setInformationForVertexPoint(EdgeSE2PointXYCov*& edge, VertexPointXYCov*& point, LaserSensorParams& params)
  {
    double beamAngle = atan2(edge->measurement()(1), edge->measurement()(0));
    ///calculate incidence angle in point frame
    Eigen::Vector2d beam = edge->measurement();

    Eigen::Vector2d normal = point->normal();
    beam = beam.normalized();
    normal = normal.normalized();
    double incidenceAngle = 0.0;
    if(point->covariance() != Eigen::Matrix2d::Identity())
      incidenceAngle = acos(normal.dot(beam));

    double d = (tan(params.angularResolution * 0.5) * edge->measurement().norm());
    ///uncertainty of the surface measurement in direction of the beam
    double dk = fabs(params.scale * (d / cos(incidenceAngle)));
    edge->information().setIdentity();
    edge->information()(0,0) = 1.0 / ((dk + params.sensorPrecision) * (dk + params.sensorPrecision));
    double cError = 0.001 * edge->measurement().norm();
    edge->information()(1,1) = 1.0 / (cError * cError);

    Eigen::Rotation2Dd rot(beamAngle);
    Eigen::Matrix2d mrot = rot.toRotationMatrix();
    edge->information() = mrot * edge->information() * mrot.transpose();
  }

}