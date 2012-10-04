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


namespace ssa{

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::CorrespondenceRejectionT()
  {
  };

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::~CorrespondenceRejectionT()
  {
  };

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  bool 
  CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(PointVertex*& reference, PointVertex*& correspondence, double& maxAngleDifference)
  {
    double angleBetweenNormals = fabs(acos(reference->globalNormal().dot(correspondence->globalNormal())));

    if(angleBetweenNormals < DEG2RAD(maxAngleDifference) )  
      return true;

    return false;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  bool 
  CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(PointVertex*& reference, PointVertex*& correspondence, double& maxAngleDifference, double& maxColorChannel)
  {
    if(correspondence->id() == reference->id())
      return false;
    if(reference->covariance() == PointMatrix::Identity() || correspondence->covariance() == PointMatrix::Identity())
      return false;

    double angleBetweenNormals = fabs(acos(reference->globalNormal().dot(correspondence->globalNormal())));
    if(angleBetweenNormals >= DEG2RAD(maxAngleDifference))
      return false;
    if(maxColorChannel == 1.0)
      return true;
    /** simple color difference check */
    double maxDiff = 255 * maxColorChannel;
    double diffCr = fabs(reference->cr - correspondence->cr);
    if(diffCr >= maxDiff)
      return false;
    double diffCg = fabs(reference->cg - correspondence->cg);
    if(diffCg >= maxDiff)
      return false;
    double diffCb = fabs(reference->cb - correspondence->cb);
    if(diffCb >= maxDiff)
      return false;
    return true;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  bool 
  CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(CorrespondenceT<PointVertex, EdgeType3>& correspondence, double& maxAngleDifference, double& maxColorChannel)
  {
    if(correspondence.cor_point->id() == correspondence.query_point->id())
      return false;
    if(correspondence.cor_point->parentVertexId() == correspondence.query_point->parentVertexId())
      return false;
    if(correspondence.query_point->covariance() == PointMatrix::Identity() || correspondence.cor_point->covariance() == PointMatrix::Identity())
      return false;

    double angleBetweenNormals = fabs(acos(correspondence.query_point->globalNormal().dot(correspondence.cor_point->globalNormal())));
    if(angleBetweenNormals >= DEG2RAD(maxAngleDifference))
      return false;
    if(maxColorChannel == 1.0)
      return true;
    /** simple color difference check */
    double maxDiff = 255 * maxColorChannel;
    double diffCr = fabs(correspondence.query_point->cr - correspondence.cor_point->cr);
    if(diffCr >= maxDiff)
      return false;
    double diffCg = fabs(correspondence.query_point->cg - correspondence.cor_point->cg);
    if(diffCg >= maxDiff)
      return false;
    double diffCb = fabs(correspondence.query_point->cb - correspondence.cor_point->cb);
    if(diffCb >= maxDiff)
      return false;
    return true;
  }
  
}