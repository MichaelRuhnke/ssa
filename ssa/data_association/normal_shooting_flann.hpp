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
  NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::NormalShootingFlann()
  {
  };

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::~NormalShootingFlann()
  {
  };

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::shootNormals(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, NormalShootingParams& params)
  {
    if(params.steps < 2){
        cerr << __PRETTY_FUNCTION__ << ": Warning params.steps is set to value < 1. Skipping normal shooting." << endl;
        return;
    }

    /** FLANN */
    std::vector<PointVertex*>& vertices = graph._verticies_points;
    PointTree kDTree;
    kDTree.copyData(vertices);
    kDTree.createKDTree();

    std::vector<int>  indices = graph.getPoseIds();
    #pragma omp parallel for schedule(dynamic, 1) shared(graph, kDTree, params)
    for(size_t i = 0; i < indices.size(); ++i)
    {
      int id = indices[i];
      shootNormals(graph, kDTree, id, params);
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::shootNormals(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::PointTree& kDTree, int& scanId, NormalShootingParams& params)
  {
    /** get vertices of scan */
    Observation<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex>& scan = graph._verticies_observations[scanId];

    #pragma omp parallel for schedule(dynamic, 20) shared(scan, graph, params)
    for(unsigned int k = 0; k < scan.size(); k=k+params.increment)
    {
      PointVertex*& point = scan[k];
      /** We apply normal shooting only if the normal is well defined */
      if(point->fixed() || point->normal().norm() == 0)
        continue;

      PointVector normal = point->normal();
      PointVertex* vp = new PointVertex;
      for(int l=-params.steps; l < params.steps; ++l)
      {
        PointVector delta;
        for(int i=0; i < Dj; ++i)
          delta(i) = ((l*params.stepSize)*normal(i));

        vp->setEstimate(point->estimate() + delta);

        //FLANN
        std::vector<int> k_indices;
        std::vector<float> k_squared_distances;

        if(kDTree.radiusSearch(vp, params.stepSize, k_indices, k_squared_distances, 1) != 0){
          PointVertex* correspondence = dynamic_cast<PointVertex*>(graph._optimizer.vertex(k_indices[0]));
          if(graph.haveCommonSensorVertex(correspondence, point))
            continue;

          if(CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(point, correspondence, params.maxAngleDifference, params.maxColorChannelDiff))
          {
            EdgeType3* edge = graph.createEdge(correspondence, point);
            #pragma omp critical
            if(edge != 0)
              graph.addEdge(edge);
          } else {
//             if(l >= 0)
//             {
//               correspondence->covariance() = PointMatrix::Identity();
//             }
          }
        }
      }
      delete vp;
    }
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::normalOutlierReject(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, NormalShootingParams& params)
  {
    if(params.steps < 2){
        cerr << __PRETTY_FUNCTION__ << ": Warning params.steps is set to value < 1. Skipping normal shooting." << endl;
        return;
    }

    /** FLANN */
    std::vector<PointVertex*>& vertices = graph._verticies_points;
    PointTree kDTree;
    kDTree.copyData(vertices);
    kDTree.createKDTree();

    std::vector<int>  indices = graph.getPoseIds();
    #pragma omp parallel for schedule(dynamic, 1) shared(graph, kDTree, params)
    for(size_t i = 0; i < indices.size(); ++i)
    {
      int id = indices[i];
      normalOutlierReject(graph, kDTree, id, params);
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::normalOutlierReject(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::PointTree& kDTree, int& scanId, NormalShootingParams& params)
  {
    /** get vertices of scan */
    Observation<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex>& scan = graph._verticies_observations[scanId];
    //FLANN
    std::vector<int> k_indices;
    std::vector<float> k_squared_distances;

    #pragma omp parallel for schedule(dynamic, 20) shared(scan, graph, params)  firstprivate(k_indices, k_squared_distances)
    for(unsigned int k = 0; k < scan.size(); ++k)
    {
      PointVertex*& point = scan[k];
      /** We apply normal shooting only if the normal is well defined */
      if(point->normal().norm() == 0.0)
        continue;

      PointVector normal = point->normal();
      PointVertex* vp = new PointVertex;
      for(int l=2; l < params.steps; ++l)
      {
        PointVector delta;
        for(int i=0; i < Dj; ++i)
          delta(i) = ((l*params.stepSize)*normal(i));

        vp->setEstimate(point->estimate() + delta);

        if(kDTree.radiusSearch(vp, params.stepSize, k_indices, k_squared_distances, 1) != 0){
          PointVertex* correspondence = dynamic_cast<PointVertex*>(graph._optimizer.vertex(k_indices[0]));
          if(graph.haveCommonSensorVertex(correspondence, point))
            continue;

          if(CorrespondenceRejectionT<EdgeType1, EdgeType2, EdgeType3>::isValid(point, correspondence, params.maxAngleDifference, params.maxColorChannelDiff))
          {
          } else {
            if(l >= 0)
            {
              correspondence->normal().setZero();
            }
          }
        }
      }
      delete vp;
    }
  }
}
