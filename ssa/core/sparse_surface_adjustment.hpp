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
  SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::SparseSurfaceAdjustmentT()
  {

  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::~SparseSurfaceAdjustmentT()
  {

  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::setGraph(SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::SSAGraph& graph)
  {
    graph_ = graph;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  typename SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::SSAGraph*
  SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::graph()
  {
    return &graph_;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::setSolver(g2o::BlockSolverX::LinearSolverType*& linearSolver)
  {
    g2o::BlockSolver< BlockSolverTraits<-1, -1> >* blockSolver = new g2o::BlockSolver< BlockSolverTraits<-1, -1> >(linearSolver);
    g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blockSolver);
    graph_._optimizer.setAlgorithm(solver);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::setVerbose(bool verbose)
  {
    graph_._optimizer.setVerbose(verbose);
    verbose_ = verbose;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::setParams(SparseSurfaceAdjustmentParams& params)
  {
    params_ = params;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  SparseSurfaceAdjustmentParams& SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::params()
  {
    return params_;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::optimizeHierarchical(int startLevel)
  {
    for(int i = startLevel; i >= 0; --i)
      optimize(i);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::optimize(int level)
  {
    double timing;
    double startTime = get_time();
    std::map<int, int> keys = graph_.getScanIds();

    for(int i = 0; i < params().ssaIterations; ++i){
      double iterationTime = get_time();
      cerr << "Iteration " << i+1 << endl;

      //! Update point covariances and normals
      timing = get_time();
      graph_.calcMeanCov(params());
      cerr << "Updating observation covariances took " << (get_time() - timing) * 1000 << " ms" << endl;

      /// Estimate data association
      timing = get_time();
        graph_.dropDataAssociation();
      /// use new abstract data association
        DataAssociationT<EdgeType1, EdgeType2, EdgeType3> da_kdtree;

	da_kdtree.setStrategy(DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::KDTREE);
        da_kdtree.apply(graph_, params(), level);

//         da_kdtree.setStrategy(DataAssociationT<EdgeType1, EdgeType2, EdgeType3>::KDTREE_REP_SUBSET);
//         da_kdtree.apply(graph_, params(), level);
//


// 	NormalShootingFlann<EdgeType1, EdgeType2, EdgeType3>::shootNormals(graph_, params().normalShooting);
        cerr << "Data association took " << (get_time() - timing) * 1000 << " ms" << endl;

      /// get edge set
      timing = get_time();
      g2o::OptimizableGraph::EdgeSet eset = graph_.getEdgesetFast();
      cerr << "Edge set construction of size " << eset.size() << " took " << (get_time() - timing) * 1000 << " ms" << endl;

      ///avoid optimizer segfault
      if(eset.size() == 0)
       continue;

      /// initializeOptimization
      timing = get_time();
      graph_._optimizer.initializeOptimization(eset);
      cerr << "initializeOptimization took " << (get_time() - timing) * 1000 << " ms" << endl;

      timing = get_time();
      if(level == 0){
	graph_._optimizer.optimize(params().g2oIterations);
      } else {
	graph_._optimizer.optimize(1);
      }
      cerr << "g2o took " << (get_time() - timing) * 1000 << " ms" << endl;

      /// moveUnoptimizedPoints
      timing = get_time();
      graph_.moveUnoptimizedPoints();
      cerr << "moveUnoptimizedPoints took " << (get_time() - timing) * 1000 << " ms" << endl;

      cerr << "Iteration took " << (get_time() - iterationTime) * 1000 << " ms" << endl;
    }
    if(params().outlierRejectionMinConnectedNeighbors > 1){
      graph()->filterOutlier(params().outlierRejectionMinConnectedNeighbors);
    }
    graph_.dropDataAssociation();
    cerr << "Computation took " << (get_time() - startTime) * 1000 << " ms" << endl;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::dumpConnectivityMatrix()
  {
    ofstream matrix("ssa-matrix-dump.txt");
      for (g2o::OptimizableGraph::EdgeSet::iterator it=graph()->_optimizer.edges().begin(); it!=graph()->_optimizer.edges().end(); it++)
      {
        matrix << (*it)->vertices()[0]->id() << " " << (*it)->vertices()[1]->id() << " 1" << endl;
      }
    matrix.close();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentT<EdgeType1, EdgeType2, EdgeType3>::optimizeColorsIncidenteWeight()
  {
    for(int i = 0; i < (int) graph()->_verticies_points.size(); ++i){
      PointVertex* point = graph()->_verticies_points[i];
      std::vector<PointVertex* > neighborhood;
      double totalWeight = 0;
      double r = 0; //point->cr;
      double g = 0; //point->cg;
      double b = 0; //point->cb;

      //find neighborhood
      for (g2o::OptimizableGraph::EdgeSet::iterator it=point->edges().begin(); it!=point->edges().end(); it++){
        EdgeType3* edge=dynamic_cast<EdgeType3*>(*it);
        if(edge){
          PointVertex* vp1=dynamic_cast<PointVertex*>(edge->vertices()[0]);
          PointVertex* vp2=dynamic_cast<PointVertex*>(edge->vertices()[1]);
          if(point->id() != vp1->id()){
            neighborhood.push_back(vp1);
          }
          if(point->id() != vp2->id()){
            neighborhood.push_back(vp2);
          }
        }
      }

      ///Add weighted color values
      for(size_t j = 0; j < neighborhood.size(); ++j){
        PointVertex*& n = neighborhood[j];
        ///find sensor edge and calculate incidence angle
        PoseVertex* nParent =  n->parentVertex();
        PointVector beam = nParent->estimate().translation() - n->estimate();
        double incidenceAngle = fabs(acos(beam.dot(n->normal())));
        double weight = 1.0; // - (incidenceAngle / (M_PI * 0.5));
//         weight = min(weight, 0.0);
//         weight *= weight;
//         if(RAD2DEG(incidenceAngle) > 60.0)
//           weight = 0.0;
        r += n->cr * weight;
        g += n->cg * weight;
        b += n->cb * weight;
        totalWeight += weight;
      }

      //set colors
      if(neighborhood.size() > 1){
        point->cr = r / totalWeight;
        point->cg = g / totalWeight;
        point->cb = b / totalWeight;
        for(size_t j = 0; j < neighborhood.size(); ++j){
          PointVertex*& n = neighborhood[j];
          n->cr = point->cr;
          n->cg = point->cg;
          n->cb = point->cb;
        }
      }
    }
  }

} //end namespace


