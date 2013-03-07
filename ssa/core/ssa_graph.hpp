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
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::SparseSurfaceAdjustmentGraphT() : maxVertexId(0), max_level_(0)
  {
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::~SparseSurfaceAdjustmentGraphT()
  {
    #pragma omp critical
    {
      clear();
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::clear(){
    _verticies_observations.clear();
    _edges_odometry.clear();
    _edges_observations.clear();
    _edges_data_association.clear();
    _edges_surface_mesh.clear();
    //std::cerr << PVAR(_optimizer.vertices().size()) << " " << PVAR(_optimizer.edges().size()) << std::endl;
    for(size_t i = 0; i < _verticies_points.size(); ++i){
      PointVertex* v = _verticies_points[i];
      _optimizer.removeVertex(v);
      _verticies_points[i] = 0;
    }
    for(size_t i = 0; i < _verticies_poses.size(); ++i){
      PoseVertex* v = _verticies_poses[i];
      _optimizer.removeVertex(v);
      _verticies_poses[i] = 0;
    }
    _verticies_poses.clear();
    _verticies_points.clear();
    //std::cerr << PVAR(_optimizer.vertices().size()) << " " << PVAR(_optimizer.edges().size()) << std::endl;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::load(std::string filename){
    std::ifstream graphstream(filename.c_str());
    if (! graphstream ){
      std::cerr << "SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::load( " << filename << ")" << ": error in loading the graph." << std::endl;
      return;
    }
    _optimizer.load(graphstream);
    std::cerr << "loaded ssa graph with " << _optimizer.vertices().size() << " vertices and " << _optimizer.edges().size() << " edges." << std::endl; 
    linkNodesToVertices();
    std::cerr << "filled access structs..." << std::endl;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::save(std::string filename){
    std::ofstream graphstream(filename.c_str());
    if (! graphstream ){
      std::cerr << "SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::save( " << filename << ")" << ": error while opening the target file." << std::endl;
      return;
    }

    g2o::Factory* factory = g2o::Factory::instance();
    for (size_t i = 0; i < _verticies_poses.size(); ++i){
      PoseVertex*& v = _verticies_poses[i];
      std::string tag = factory->tag(v);
      if (tag.size() > 0) {
        graphstream << tag << " " << v->id() << " ";
        v->write(graphstream);
        graphstream << std::endl;
      }
    }

    for (size_t i = 0; i < _verticies_points.size(); ++i){
      PointVertex*& v = _verticies_points[i];
      std::string tag = factory->tag(v);
      if (tag.size() > 0) {
        graphstream << tag << " " << v->id() << " ";
        v->write(graphstream);
        graphstream << std::endl;
      }
    }

    for (size_t i = 0; i < _edges_odometry.size(); ++i){
      SLAMEdgeType*& e = _edges_odometry[i];
      std::string tag = factory->tag(e);
      if (tag.size() > 0) {
        graphstream << tag << " ";
        if (e->id() >= 0)
          graphstream << e->id() << " ";
        for (std::vector<HyperGraph::Vertex*>::const_iterator it = e->vertices().begin(); it != e->vertices().end(); ++it) {
          OptimizableGraph::Vertex* v = static_cast<OptimizableGraph::Vertex*>(*it);
          graphstream << v->id() << " ";
        }
        //graphstream << e->level() << " ";
        e->write(graphstream);
        graphstream << std::endl;
      }
    }

    for (size_t i = 0; i < _edges_observations.size(); ++i){
      SensorEdgeType*& e = _edges_observations[i];
      std::string tag = factory->tag(e);
      if (tag.size() > 0) {
        graphstream << tag << " ";
        if (e->id() >= 0)
          graphstream << e->id() << " ";
        for (std::vector<HyperGraph::Vertex*>::const_iterator it = e->vertices().begin(); it != e->vertices().end(); ++it) {
          OptimizableGraph::Vertex* v = static_cast<OptimizableGraph::Vertex*>(*it);
          graphstream << v->id() << " ";
        }
        //graphstream << e->level() << " ";
        e->write(graphstream);
        graphstream << std::endl;
      }
    }

    for (size_t i = 0; i < _edges_data_association.size(); ++i){
      DataAssociationEdgeType*& e = _edges_data_association[i];
      std::string tag = factory->tag(e);
      if (tag.size() > 0) {
        graphstream << tag << " ";
        if (e->id() >= 0)
          graphstream << e->id() << " ";
        for (std::vector<HyperGraph::Vertex*>::const_iterator it = e->vertices().begin(); it != e->vertices().end(); ++it) {
          OptimizableGraph::Vertex* v = static_cast<OptimizableGraph::Vertex*>(*it);
          graphstream << v->id() << " ";
        }
        //graphstream << e->level() << " ";
        e->write(graphstream);
        graphstream << std::endl;
      }
    }

//     set<Vertex*, VertexIDCompare> verticesToSave;
//     for (HyperGraph::EdgeSet::const_iterator it = _optimizer.edges().begin(); it != _optimizer.edges().end(); ++it) {
//       OptimizableGraph::Edge* e = static_cast<OptimizableGraph::Edge*>(*it);
//       for (vector<HyperGraph::Vertex*>::const_iterator it = e->vertices().begin(); it != e->vertices().end(); ++it) {
//         verticesToSave.insert(static_cast<OptimizableGraph::Vertex*>(*it));
//       }
//     }
// 
//     for (set<Vertex*, VertexIDCompare>::const_iterator it = verticesToSave.begin(); it != verticesToSave.end(); ++it){
//       OptimizableGraph::Vertex* v = *it;
//       _optimizer.saveVertex(graphstream, v);
//     }
// 
//     EdgeContainer edgesToSave;
//     for (HyperGraph::EdgeSet::const_iterator it = _optimizer.edges().begin(); it != _optimizer.edges().end(); ++it) {
//       const OptimizableGraph::Edge* e = dynamic_cast<const OptimizableGraph::Edge*>(*it);
//       if (e->level() == level)
//         edgesToSave.push_back(const_cast<Edge*>(e));
//     }
//     sort(edgesToSave.begin(), edgesToSave.end(), EdgeIDCompare());
// 
//     for (EdgeContainer::const_iterator it = edgesToSave.begin(); it != edgesToSave.end(); ++it) {
//       OptimizableGraph::Edge* e = *it;
//       _optimizer.saveEdge(graphstream, e);
//     }

    //_optimizer.save(graphstream);
    graphstream.close();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::linkNodesToVertices(){

  _verticies_poses.clear();
  _verticies_points.clear();
  _verticies_observations.clear();
  _edges_odometry.clear();
  _edges_observations.clear();
  _edges_data_association.clear();
  _edges_surface_mesh.clear();

    for (g2o::OptimizableGraph::VertexIDMap::const_iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
      PoseVertex* v=dynamic_cast<PoseVertex*>(it->second);
      if(v){  
          //Fix first pose vertex
          if(_verticies_poses.size() == 0)
            v->setFixed(true);
          _verticies_poses.push_back(v);
          maxVertexId = max(maxVertexId, v->id());
      }
    }

    for (g2o::OptimizableGraph::VertexIDMap::const_iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
      PointVertex* v=dynamic_cast<PointVertex*>(it->second);
      if(v){  
        maxVertexId = max(maxVertexId, v->id());
        v->setParentVertex(dynamic_cast<PoseVertex*> (_optimizer.vertices()[v->parentVertexId()]));
        if(v->parentVertex()){
          _verticies_points.push_back(v);
         }
      }
    }

    for (g2o::OptimizableGraph::EdgeSet::iterator it=_optimizer.edges().begin(); it!=_optimizer.edges().end(); it++){
      SLAMEdgeType* e1=dynamic_cast<SLAMEdgeType*>(*it);
      double chi2 = 0.0;
      if(e1){
        e1->computeError();
        chi2 = e1->chi2();
        if(chi2 < 0.0){
          std::cerr << "Error: NEGATIVE INFORMATION MATRIX... THIS SHOULD NOT HAPPEN..." << std::endl;
          std::cerr << "SLAMEdgeType " << chi2 << "\t";
          std::cerr << e1->information() << std::endl;
          e1->information() = -1.0 * e1->information();
          //exit(-1);
        }
        _edges_odometry.push_back(e1);
      } else {
        SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
        if(e2){
          e2->computeError();
          chi2 = e2->chi2();
          if(chi2 < 0.0){
            std::cerr << "Error: NEGATIVE INFORMATION MATRIX... THIS SHOULD NOT HAPPEN..." << " SensorEdgeType " << e2->vertices()[0]->id() << " " << e2->vertices()[1]->id()  << " " << chi2 << std::endl;
            e2->information().setIdentity();
          }
          _edges_observations.push_back(e2);
          max_level_ = max(max_level_, e2->level());
        } else {
          DataAssociationEdgeType* e3=dynamic_cast<DataAssociationEdgeType*>(*it);
          if(e3){
            e3->computeError();
            chi2 = e3->chi2();
            if(chi2 < 0.0){
              std::cerr << "Error: NEGATIVE INFORMATION MATRIX... THIS SHOULD NOT HAPPEN... " << e3->vertices()[0]->id() << " " << e3->vertices()[1]->id() << " " << std::endl;
//               std::cerr << std::fixed << chi2 << "\r\n";
//               std::cerr << e3->information() << std::endl;
              _optimizer.removeEdge(e3);
//               exit(-1);
            } else {
              _edges_surface_mesh.push_back(e3);
            }
          }
        }
      }
    }
    fillObservations();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::fillObservations()
  {
    _verticies_observations.clear();
    _edges_observations.clear();
    //typename
    for (typename std::vector<PoseVertex * >::iterator it=_verticies_poses.begin(); it!=_verticies_poses.end(); it++){
      PoseVertex*& v=(*it);
      if(v){
        for (g2o::OptimizableGraph::EdgeSet::iterator itt=v->edges().begin(); itt!=v->edges().end(); itt++){
          SensorEdgeType* edge=dynamic_cast<SensorEdgeType*>(*itt);
          if(edge){
            PoseVertex* v1=dynamic_cast<PoseVertex*>(edge->vertices()[0]);
            PointVertex* v2=dynamic_cast<PointVertex*>(edge->vertices()[1]);
            if(v1 && v2){
              _verticies_observations[v->id()].push_back(v2);
              //v2->neighbors.clear();
            }
            _edges_observations.push_back(edge);
          }
        }
      }
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  std::map<int, int> SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getScanIds(){
    std::cerr << __PRETTY_FUNCTION__ << ": deprecated and slow use getPoseIds() instead." << std::endl;
    std::map<int, int> keys;
    for (g2o::OptimizableGraph::VertexIDMap::const_iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
      PoseVertex* v=dynamic_cast<PoseVertex*>(it->second);
      if(v)
        keys[v->id()] = v->id();
    }
    return keys;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  std::vector<int> SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getPoseIds(){
    std::vector<int>  indices;
    //for (g2o::OptimizableGraph::VertexIDMap::const_iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
    indices.reserve(_verticies_poses.size());
    for(size_t i = 0; i < _verticies_poses.size(); ++i){
      indices.push_back(_verticies_poses[i]->id());
    }
    return indices;
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::calcMeanCov(SparseSurfaceAdjustmentParams& params){
    std::vector<int> keys = getPoseIds();
    #pragma omp parallel for shared(params, keys) 
    for(int i = 0; i < (int) keys.size(); ++i)
      calcMeanCov(params, keys[i]);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::calcMeanCov(SparseSurfaceAdjustmentParams& params, int scanId){
    Observation<PointVertex>& observation = getObservationOfScan(scanId);
    for(size_t j = 0; j < observation.size(); ++j){
      neighborCache().getCovariance(observation[j]->id(), observation[j]->covariance(), params);
      observation[j]->updateNormal(observation[j]->covariance());
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::fillNeighborCache(SparseSurfaceAdjustmentParams& params)
  {
    neighborCache().clear();
    std::vector<int> keys = getPoseIds();
    //#pragma omp parallel for shared(params, keys) num_threads(params.maxThreads)
    for(size_t i = 0; i < keys.size(); ++i){
      fillNeighborCache(params, keys[i]);
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::fillNeighborCache(SparseSurfaceAdjustmentParams& params, int scanId)
  {
    Observation<PointVertex>& observation = getObservationOfScan(scanId);
    if(observation.size() > 0)
      neighborCache().createFromScan(observation, params);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  g2o::OptimizableGraph::EdgeSet SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getEdgeset()
  {
    g2o::OptimizableGraph::EdgeSet eset;
    std::map<int, bool>  pointWithCorrespondence;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=_optimizer.edges().begin(); it!=_optimizer.edges().end(); it++){
        SLAMEdgeType* e1=dynamic_cast<SLAMEdgeType*>(*it);
        if(e1){
           eset.insert(e1);
        } else {
          DataAssociationEdgeType* e3=dynamic_cast<DataAssociationEdgeType*>(*it);
          if(e3){
            PointVertex* v1=dynamic_cast<PointVertex*>(e3->vertices()[0]);
            PointVertex* v2=dynamic_cast<PointVertex*>(e3->vertices()[1]);
            pointWithCorrespondence[v1->id()] = true;
            pointWithCorrespondence[v2->id()] = true;
            eset.insert(e3);
          }
        }
    }
    for (g2o::OptimizableGraph::EdgeSet::iterator it=_optimizer.edges().begin(); it!=_optimizer.edges().end(); it++){
      SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
      if(e2){
        PoseVertex* v=dynamic_cast<PoseVertex*>(e2->vertices()[0]);
        PointVertex* v1=dynamic_cast<PointVertex*>(e2->vertices()[1]);
        if(pointWithCorrespondence[v1->id()]){
          eset.insert(e2);
        } else{
          if(!v->fixed()){
            _edges_points_to_move.push_back(e2);
            PointVector backupTransform = v->estimate().inverse() * v1->estimate();
            _edges_points_to_move_transforms.push_back(backupTransform);
          }
        }
      } 
    }
    return eset;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  g2o::OptimizableGraph::EdgeSet SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getEdgesetFast()
  {
    g2o::OptimizableGraph::EdgeSet eset;
    for (size_t i = 0; i < _edges_odometry.size(); ++i){
      eset.insert(_edges_odometry[i]);
    }
    std::map<int, bool> isAssigned;
    for (size_t i = 0; i < _edges_data_association.size(); ++i){
      DataAssociationEdgeType*& e3=_edges_data_association[i];
      isAssigned[e3->vertices()[0]->id()] = true;
      isAssigned[e3->vertices()[1]->id()] = true;
      eset.insert(e3);
    }

//     #pragma omp parallel for schedule(dynamic, 2) shared(eset) 
    for (size_t i = 0; i < _edges_observations.size(); ++i){
      SensorEdgeType*& e2=_edges_observations[i];
      PoseVertex* v=static_cast<PoseVertex*>(e2->vertices()[0]);
      PointVertex* v1=static_cast<PointVertex*>(e2->vertices()[1]);
        if(isAssigned[v1->id()]){
//         #pragma omp critical
          eset.insert(e2);
        } else {
          PointVector backupTransform = v->estimate().inverse() * v1->estimate();
//           #pragma omp critical
          {
            _edges_points_to_move.push_back(e2);
            _edges_points_to_move_transforms.push_back(backupTransform);
          }
        }
    }

    for (size_t i = 0; i < _edges_data_association_gicp.size(); ++i){
      eset.insert(_edges_data_association_gicp[i]);
    }
//     std::cerr << "SSA p-t-p correspondences       \t " << _edges_data_association.size() << std::endl;
//     std::cerr << "GICP p-t-p correspondences      \t " << _edges_data_association_gicp.size() << std::endl;
    std::cerr << "Points keep out of optimization \t " << _edges_points_to_move.size() << std::endl;
    std::cerr << "Number of edges for optimization\t " << eset.size() << std::endl;
//     std::cerr << PVAR(_edges_points_to_move.size()) << " " << PVAR(eset.size()) << std::endl;
    return eset;
  }


//   template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
//   g2o::OptimizableGraph::EdgeSet SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getEdgesetMesh()
//   {
//     g2o::OptimizableGraph::EdgeSet eset;
//     std::map<int, bool>  pointWithCorrespondence;
//     //TODO: switch to _edges_surface_mesh instead (should be faster)
//     for (g2o::OptimizableGraph::EdgeSet::iterator it=_optimizer.edges().begin(); it!=_optimizer.edges().end(); it++){
//         DataAssociationEdgeType* e3=dynamic_cast<DataAssociationEdgeType*>(*it);
//         if(e3){
//           PointVertex* v1=dynamic_cast<PointVertex*>(e3->vertices()[0]);
//           PointVertex* v2=dynamic_cast<PointVertex*>(e3->vertices()[1]);
//           //TODO: if(e3->isMeshEdge)
//           if(e3->isMeshEdge){
//             pointWithCorrespondence[v1->id()] = true;
//             pointWithCorrespondence[v2->id()] = true;
//             eset.insert(e3);
//           } 
//         }
//     }
//     for (g2o::OptimizableGraph::EdgeSet::iterator it=_optimizer.edges().begin(); it!=_optimizer.edges().end(); it++){
//       SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
//       if(e2){
//         PointVertex* v1=dynamic_cast<PointVertex*>(e2->vertices()[1]);
//         if(pointWithCorrespondence[v1->id()]){
//           eset.insert(e2);
//         }
//       } 
//     }
//     return eset;
//   }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::moveUnoptimizedPoints(){
    #pragma omp parallel for schedule(dynamic, 20) 
    for(size_t i = 0; i < _edges_points_to_move.size(); ++i){
      SensorEdgeType*& e=_edges_points_to_move[i];
      PoseVertex*    v1=static_cast<PoseVertex*>(e->vertices()[0]);
      PointVertex* v2=static_cast<PointVertex*>(e->vertices()[1]);
      v2->setEstimate(v1->estimate() * _edges_points_to_move_transforms[i]);
    }
    //std::cerr << "moved " << _edges_points_to_move.size() << " point without correspondences." << std::endl;
    _edges_points_to_move.clear();
    _edges_points_to_move_transforms.clear();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::filterOutlier(int minNeighbors=2){
    for(size_t i = 0; i < _verticies_points.size(); ++i){
      PointVertex* v =_verticies_points[i];
      int neighbors = v->edges().size();
//       for (g2o::OptimizableGraph::EdgeSet::iterator it=v->edges().begin(); it!=v->edges().end(); it++){
//         DataAssociationEdgeType* edge=dynamic_cast<DataAssociationEdgeType*>(*it);
//         if(edge != 0)
//           neighbors++;
//       }
      if(neighbors < minNeighbors){
          v->covariance() = PointMatrix::Identity();
      } else {
        //cerr << "neighbors " << neighbors << std::endl;
      }
    }
  }

//   template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
//   typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointTree 
//   SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getKDTreeOfScans
//   (int& scanId)
//   {
//     PointTree kdTree;
//     std::map<int, int> keys = getScanIds();
//     for(std::map<int, int>::iterator it = keys.begin(); it != keys.end(); ++it){
//       if(it->second != scanId){
// 
//         Observation<PointVertex>& o = _verticies_observations[it->second];
//         for(unsigned int j = 0; j < o.size(); ++j){
//           PointVertex* point = o[j];
//           kdTree.insert(point);
//         }
//       }
//     }
//     kdTree.optimize();
//     return kdTree;
//   }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  std::vector<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex* >
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getPointVerticesOfScans
  (int& scanId)
  {
    std::vector<PointVertex* > vertices;
    vertices.reserve(_verticies_points.size());

    std::vector<int> keys = getPoseIds();
    for(size_t i = 0; i < keys.size(); ++i){
      if(keys[i] != scanId){
        Observation<PointVertex>& o = _verticies_observations[keys[i]];
        for(unsigned int j = 0; j < o.size(); ++j){
          PointVertex* point = o[j];
          vertices.push_back(point);
        }
      }
    }
    return vertices;
  }

//   template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
//   typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointTree 
//   SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getKDTreeOfScan
//   (int& scanId)
//   {
//     PointTree kdTree;
//     Observation<PointVertex>& o = _verticies_observations[scanId];
//     for(unsigned int j = 0; j < o.size(); ++j){
//       PointVertex* point = o[j];
//       kdTree.insert(point);
//     }
//     kdTree.optimize();
//     return kdTree;
//   }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  std::vector<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex* >
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getPointVerticesOfScan
  (int& scanId)
  {
    std::vector<PointVertex* > vertices;
    Observation<PointVertex>& o = _verticies_observations[scanId];
    for(unsigned int j = 0; j < o.size(); ++j){
      PointVertex* point = o[j];
      vertices.push_back(point);
    }
    return vertices;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  Observation<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex>&
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getObservationOfScan
  (int& scanId)
  {
    return _verticies_observations[scanId];
  }
  
  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  int SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addVertex(PoseVertex*& v)
  {
    /** Handling for automatic vertex ids */
    maxVertexId = std::max(maxVertexId, v->id());
    if(maxVertexId == 0)
      v->setFixed(true);
    if(v->id() == 0) /** if vertex id is not set */
      v->setId(maxVertexId++);
    /** add to optimizer and internal pose struct */
    _optimizer.addVertex(v);
    _verticies_poses.push_back(v);
    return v->id();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  int SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addVertex(PointVertex*& v)
  {
    /** Handling for automatic vertex ids */
    maxVertexId = std::max(maxVertexId, v->id());
    if(v->id() == 0) /** if vertex id is not set */
      v->setId(maxVertexId++);
    /** add to optimizer and internal pose struct */
    _optimizer.addVertex(v);
    _verticies_points.push_back(v);
    _verticies_observations[v->parentVertex()->id()].push_back(v);

    /** force recomputation of index 
        TODO: this could be implemented in a more elegant way! */
    _hierarchical_point_vertex_index.clear();
    return v->id();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  int SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addVertex(PointVertex*& v, std::vector< PointVertex* >& neighbors)
  {
    int id = addVertex(v);
    cached_neighbors_.addVertex(v, neighbors);
    return v->id();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::removeVertex(PoseVertex*& v)
  {
    for (g2o::OptimizableGraph::EdgeSet::iterator it=v->edges().begin(); it!=v->edges().end(); it++){
      SLAMEdgeType* e1=dynamic_cast<SLAMEdgeType*>(*it);
      if(e1){
        removeEdge(e1);
      } else {
        SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
        if(e2){
          removeEdge(e2);
        }
      }
    }
    typename std::vector<PoseVertex* >::iterator it = std::find(_verticies_poses.begin(), _verticies_poses.end(), v);
    _verticies_poses.erase(it);
    _optimizer.removeVertex(v);

  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::removeVertex(PointVertex*& v)
  {

    for (g2o::OptimizableGraph::EdgeSet::iterator it=v->edges().begin(); it!=v->edges().end(); it++){
      SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
      if(e2){
        removeEdge(e2);
      } else {
        DataAssociationEdgeType* e3=dynamic_cast<DataAssociationEdgeType*>(*it);
        if(e3){
          removeEdge(e3);
        }
      }
    }

    std::vector<PoseVertex* > parents = v->parentVertices();
    for(typename std::vector<PoseVertex* >::const_iterator it = parents.begin(); it != parents.end(); ++it)
    {
      Observation<PointVertex>& observation = _verticies_observations[(*it)->id()];
      observation.removeVertex(v);
    }

    typename std::vector<PointVertex* >::iterator itt = std::find(_verticies_points.begin(), _verticies_points.end(), v);
    if(itt != _verticies_points.end())
      _verticies_points.erase(itt);

    _optimizer.removeVertex(v);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addEdge(SLAMEdgeType*& edge)
  {
    _edges_odometry.push_back(edge);
    _optimizer.addEdge(edge);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addEdge(SensorEdgeType*& edge)
  {
    _edges_observations.push_back(edge);
    _optimizer.addEdge(edge);
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addEdge(DataAssociationEdgeType*& edge)
  {
    if(edge == 0)
      return;
    _edges_data_association.push_back(edge);
    _optimizer.addEdge(edge);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addEdge(g2o::Edge_V_V_GICP*& edge)
  {
//     std::cerr << __PRETTY_FUNCTION__ << std::endl;
    if(edge == 0)
      return;

    _edges_data_association_gicp.push_back(edge);
//     std::cerr << _edges_data_association_gicp.size() << std::endl;
    _optimizer.addEdge(edge);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::addEdge(PointVertex*& from, PointVertex*& to)
  {
    if(from == to){
      std::cerr << "Warning: you tried to create an edge between an observation and itself." << std::endl;
      return;
    }

    for (g2o::OptimizableGraph::EdgeSet::iterator it=from->edges().begin(); it!=from->edges().end(); it++){
      DataAssociationEdgeType* e1=dynamic_cast<DataAssociationEdgeType*>(*it);
      if(e1){
        PointVertex* v1=dynamic_cast<PointVertex*>(e1->vertices()[0]);
        PointVertex* v2=dynamic_cast<PointVertex*>(e1->vertices()[1]);
        if((v1->id() == from->id() && v2->id() == to->id()) || (v1->id() == to->id() && v2->id() == from->id())){
          //cerr << "Edge exists already..." << std::endl;
          return;
        }
      }
    }

    DataAssociationEdgeType* e = new DataAssociationEdgeType;
    PointVector v;
    v.fill(0.0);
    e->setMeasurement(v);
    e->vertices()[0]=from;
    e->vertices()[1]=to;
    e->information() = PointMatrix::Identity();

    e->computeError();
    addEdge(e);
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::DataAssociationEdgeType* SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::createEdge(PointVertex*& from, PointVertex*& to)
  {
    if(from == to){
      std::cerr << "Warning: you tried to create an edge between an observation and itself." << std::endl;
      return 0;
    }

    for (g2o::OptimizableGraph::EdgeSet::iterator it=from->edges().begin(); it!=from->edges().end(); it++){
      DataAssociationEdgeType* e1=dynamic_cast<DataAssociationEdgeType*>(*it);
      if(e1){
        PointVertex* v1=dynamic_cast<PointVertex*>(e1->vertices()[0]);
        PointVertex* v2=dynamic_cast<PointVertex*>(e1->vertices()[1]);
        if((v1->id() == from->id() && v2->id() == to->id()) || (v1->id() == to->id() && v2->id() == from->id())){
          //cerr << "Edge exists already..." << std::endl;
          return 0;
        }
      }
    }

    DataAssociationEdgeType* e = new DataAssociationEdgeType;
    PointVector v = (to->estimate() - from->estimate()) * 0.5;
    v.fill(0.0);
    e->setMeasurement(v);
    e->vertices()[0]=from;
    e->vertices()[1]=to;
    e->information() = PointMatrix::Identity();
    e->computeError();
    return e;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::removeEdge(SLAMEdgeType* edge)
  {
    if(edge == 0)
      return;

   typename std::vector<SLAMEdgeType* >::iterator it = std::find(_edges_odometry.begin(), _edges_odometry.end(), edge);
   if(it != _edges_odometry.end())
     _edges_odometry.erase(it);

   _optimizer.removeEdge(edge);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::removeEdge(SensorEdgeType* edge)
  {
    if(edge == 0)
      return;

   typename std::vector<SensorEdgeType* >::iterator it = std::find(_edges_observations.begin(), _edges_observations.end(), edge);
   if(it != _edges_observations.end())
     _edges_observations.erase(it);

   _optimizer.removeEdge(edge);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::removeEdge(DataAssociationEdgeType* edge)
  {
    if(edge == 0)
      return;

   typename std::vector<DataAssociationEdgeType* >::iterator it = std::find(_edges_data_association.begin(), _edges_data_association.end(), edge);
   if(it != _edges_data_association.end())
     _edges_data_association.erase(it);

   _optimizer.removeEdge(edge);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::dropDataAssociation()
  {

    for(size_t i = 0; i < _edges_data_association.size(); ++i){
      DataAssociationEdgeType*& e = _edges_data_association[i];
      if(e){
        _optimizer.removeEdge(e);
      }
    }
    _edges_data_association.clear();

    for(size_t i = 0; i < _edges_data_association_gicp.size(); ++i){
      g2o::Edge_V_V_GICP*& e = _edges_data_association_gicp[i];
      if(e){
        _optimizer.removeEdge(e);
      }
    }
    _edges_data_association_gicp.clear();
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::dropSensorEdges()
  {
    int j = 0;
    for(size_t i = 0; i < _edges_observations.size(); ++i){
      SensorEdgeType* e = _edges_observations[i];
      _edges_observations[i] = 0;
      if(e){
        _optimizer.removeEdge(e);
        j++;
      }
    }
    //std::cerr << "deleted " << j << " edges. " << _edges_observations.size() << std::endl;
    _edges_observations.clear();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::dropSLAMEdges()
  {

    for(size_t i = 0; i < _edges_odometry.size(); ++i){
      SLAMEdgeType*& e = _edges_odometry[i];
      if(e){
        _optimizer.removeEdge(e);
      }
    }
    _edges_odometry.clear();
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::unfixVertices()
  {
    setFixedVertices(false);
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::setFixedVertices(bool fixed)
  {
    bool fixedOneVertex = false;
    for (g2o::OptimizableGraph::VertexIDMap::iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
      PoseVertex* p = dynamic_cast<PoseVertex* >(it->second);
      if(p){
        if(!fixedOneVertex){ //the first pose vertex has to be fixed 
          p->setFixed(true);
          fixedOneVertex = true;
        } else {
          p->setFixed(fixed);
        }
      } else {
        PointVertex* vp=dynamic_cast<PointVertex*>(it->second);
        if(vp){
          vp->setFixed(fixed);
        }
      } 
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::fixDataAssociation()
  {

//     for(size_t i = 0; i < _edges_data_association.size(); ++i){
//       DataAssociationEdgeType*& e = _edges_data_association[i];
//       if(e){
//         e->isMeshEdge = true;
//         _edges_surface_mesh.push_back(e);
//       }
//     }
//     _edges_data_association.clear();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::pruneGraph(double minDistance)
  {
    std::cerr << PVAR(_optimizer.vertices().size()) << std::endl;
    std::map<int, bool> deletedVertex;
    std::vector<PointVertex* > verticesToDelete;
    double timing = get_time();

    std::map<double, DataAssociationEdgeType*> da_edges;
    int count = 0;
    int count2 = 0;
    for(size_t i = 0; i < _edges_data_association.size(); ++i){
      DataAssociationEdgeType*& e = _edges_data_association[i];
      if(e){
        PointVertex* v1=static_cast<PointVertex*>(e->vertices()[0]);
        PointVertex* v2=static_cast<PointVertex*>(e->vertices()[1]);
        if(!v1 || !v2 || deletedVertex[v2->id()] || deletedVertex[v1->id()] || v1->id() == v2->id())
            continue;
        double distance = (v1->estimate() - v2->estimate()).norm();
        if(distance <= minDistance && !haveCommonSensorVertex(v1, v2)){
          da_edges[distance + 1e-10*v1->id()] = e;  
          count++;
          _fixedDataAssociation[v1->parentVertex()->id()][v2->parentVertex()->id()] = true;
          _fixedDataAssociation[v2->parentVertex()->id()][v1->parentVertex()->id()] = true;
        } else {
          count2++;
        }
      }
    }
    std::cerr << PVAR(count) << " " << PVAR(count2) << " " << PVAR(da_edges.size()) << std::endl;

    for(typename std::map<double, DataAssociationEdgeType*>::iterator it = da_edges.begin(); it != da_edges.end(); ++it){
      DataAssociationEdgeType*& e = it->second;
      if(!e)
        continue;

      PointVertex* v1=static_cast<PointVertex*>(e->vertices()[0]);
      PointVertex* v2=static_cast<PointVertex*>(e->vertices()[1]);

      if(v1->edges().size() < v2->edges().size())
        std::swap(v1,v2);

      for (g2o::HyperGraph::EdgeSet::iterator it=v2->edges().begin(); it!=v2->edges().end(); it++){
        SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
        if(e2){
          /** Update Observation */
          Observation<PointVertex>& observation = _verticies_observations[e2->vertices()[0]->id()];
          observation.push_back(v1);
          observation.removeVertex(v2);
          /** Update Edge */
          e2->vertices()[1] = v1;
          v1->edges().insert(e2);
          v2->edges().erase(it);
        }
      }
      /** Mark Vertex as deleted */
      deletedVertex[v2->id()] = true;
      verticesToDelete.push_back(v2);
    }
   std::cerr << "iterating through DA edges took " << (get_time() - timing) * 1000 << " ms." << std::endl;


   timing = get_time();
   dropDataAssociation();
   std::cerr << "Drop DA edges took " << (get_time() - timing) * 1000 << " ms." << std::endl;

  std::cerr << "pruned " << verticesToDelete.size() << " vertices." << std::endl;
//   timing = get_time();
//   for(size_t i = 0; i < _verticies_points.size(); ++i){
//     PointVertex* v = _verticies_points[i];
//     if(v->covariance() == PointMatrix::Identity()){
//       verticesToDelete.push_back(v);
//     }
//   }
//   std::cerr << "pruned " << verticesToDelete.size() << " vertices." << std::endl;
  std::cerr << "Removing outlier took " << (get_time() - timing) * 1000 << " ms." << std::endl;
  timing = get_time();
//   for(size_t i = 0; i < verticesToDelete.size(); ++i){
//     PointVertex* p = verticesToDelete[i];
//     if(p)
//       removeVertex(p);
//     verticesToDelete[i] = 0;
//   }

  std::cerr << PVAR(_optimizer.vertices().size()) << std::endl;
  fillObservations();
  std::cerr << "deleting vertices took " << (get_time() - timing) * 1000 << " ms." << std::endl;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::marginalizeEdges(PoseVertex* v)
  {
/*  ///EdgeLabeler not yet released
    EdgeLabeler labeler(&_optimizer);
    //std::set<SLAMEdgeType*> slamEdges;
    std::set<OptimizableGraph::Edge*> slamEdges;
    for (g2o::OptimizableGraph::VertexIDMap::iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
      PoseVertex* p = dynamic_cast<PoseVertex* >(it->second);
      if(p){
        if(p->id() == v->id())
          continue;
        if(!connected(p, v)){
          //cerr << "not connected " << p->id() << " and " << v->id() << std::endl;
          continue;
        }
        SLAMEdgeType*  edge = new SLAMEdgeType;
        edge->vertices()[0]=p;
        edge->vertices()[1]=v;
        addEdge(edge);
        slamEdges.insert(edge);
      } 
    }
    int error = labeler.labelEdges(slamEdges);
    if(error == -1)
        std::cerr << "Warning: labelEdges failed..." << std::endl;

//     for (g2o::OptimizableGraph::VertexIDMap::iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
//       PoseVertex* p = dynamic_cast<PoseVertex* >(it->second);
//       if(p){
//         p->setFixed(true);
//       } else {
//         PointVertex* vp=dynamic_cast<PointVertex*>(it->second);
//         if(vp){
//           vp->setFixed(true);
//         }
//       } #
//     }
    dropDataAssociation(); */
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::pruneDataAssociation(int maxEdges)
  {

    int count = 0;
//     std::vector<DataAssociationEdgeType* > edgesToRemove;
    _edges_data_association.clear();
    for(size_t i = 0; i < _verticies_points.size(); ++i){
      PointVertex* vp = _verticies_points[i];
      if(vp){
        if(vp->edges().size() <= 1)
          continue;
        typename std::deque<DataAssociationEdgeType*> edgesOfVertex;
        for (g2o::OptimizableGraph::EdgeSet::iterator itt=vp->edges().begin(); itt!=vp->edges().end(); itt++){
          DataAssociationEdgeType* e2=dynamic_cast<DataAssociationEdgeType*>(*itt);
          if(e2)
            edgesOfVertex.push_back(e2);
        }
        std::sort(edgesOfVertex.begin(), edgesOfVertex.end(), DataAssociationEdgeType::DataAssociationEdgeComp);
        count = 0;
        for(typename std::deque<DataAssociationEdgeType*>::const_iterator it=edgesOfVertex.begin(); it!=edgesOfVertex.end(); ++it){
          DataAssociationEdgeType* e2 = (*it);
          if(e2){
            if(count >= maxEdges && (int) e2->vertices()[0]->edges().size() > (maxEdges + 1) && (int) e2->vertices()[1]->edges().size() > (maxEdges + 1)){
              _optimizer.removeEdge(e2);
              count++;
            } else {
              _edges_data_association.push_back(e2);
            }
            count++;
          }
        }


      }
      std::cerr << " " << (double) i / (double) _verticies_points.size() << "."<< std::endl;
    }
//     for(size_t i = 0; i < edgesToRemove.size(); ++i){
//       DataAssociationEdgeType* e2 = edgesToRemove[i];
//       _optimizer.removeEdge(e2);
//     }
    std::cerr << " done! " << std::endl;
    linkNodesToVertices();
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  bool SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::connected(PoseVertex*& v1, PoseVertex*& v2)
  {
    for (g2o::OptimizableGraph::EdgeSet::iterator it=v1->edges().begin(); it!=v1->edges().end(); it++){
      SensorEdgeType* e=dynamic_cast<SensorEdgeType*>(*it);
      if(e){
        PointVertex* vp = dynamic_cast<PointVertex* >(e->vertices()[1]);
        if(vp){
          for (g2o::OptimizableGraph::EdgeSet::iterator itt=vp->edges().begin(); itt!=vp->edges().end(); itt++){
            DataAssociationEdgeType* e2=dynamic_cast<DataAssociationEdgeType*>(*itt);
            if(e2){
              PointVertex* vp2 = dynamic_cast<PointVertex* >(e2->vertices()[0]);
              PointVertex* vp3 = dynamic_cast<PointVertex* >(e2->vertices()[1]);
              //check direction 
              if(vp2 && vp3)
                if(vp2->id() == vp->id())
                  vp2 = vp3;
            
              if(vp2){
                for (g2o::OptimizableGraph::EdgeSet::iterator ittt=vp2->edges().begin(); ittt!=vp2->edges().end(); ittt++){
                  SensorEdgeType* e3=dynamic_cast<SensorEdgeType*>(*ittt);
                  if(e3){
                    PoseVertex* p = dynamic_cast<PoseVertex* >(e3->vertices()[0]);
                    if(p)
                      if(v2->id() == p->id())
                        return true;
                  }
                }
              }
            }
          }
        }
      }
    }
    return false;
  }


  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  bool SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::connected(PointVertex*& v1, PointVertex*& v2)
  {
    for (g2o::OptimizableGraph::EdgeSet::iterator it=v1->edges().begin(); it!=v1->edges().end(); it++){
      DataAssociationEdgeType* e2=dynamic_cast<DataAssociationEdgeType*>(*it);
      if(e2){
        if((e2->vertices()[0] == v1 && e2->vertices()[1] == v2)|| (e2->vertices()[0] == v2 && e2->vertices()[1] == v1))
          return true;
      }
    }
    return false;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  bool SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::haveCommonSensorVertex(PointVertex*& v1, PointVertex*& v2)
  {
    if(v1->id()  == v2->id())
      return true;
    bool haveCommonPose = false;
    std::map< int, bool> parentMatch;

    for (g2o::OptimizableGraph::EdgeSet::iterator it=v1->edges().begin(); it!=v1->edges().end(); it++){
      SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
      if(e2){
        if(!e2->vertices()[0]){
          std::cerr << __PRETTY_FUNCTION__ << " vertex not available. This SHOULD not happen! " << std::endl;
          exit(0);
        }
        parentMatch[e2->vertices()[0]->id()] = true;
      }
    }

    for (g2o::OptimizableGraph::EdgeSet::iterator it=v2->edges().begin(); it!=v2->edges().end(); it++){
      SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
      if(e2){
        if(!e2->vertices()[0]){
          std::cerr << __PRETTY_FUNCTION__ << " vertex not available. This SHOULD not happen! " << std::endl;
          exit(0);
        }
        if(parentMatch[e2->vertices()[0]->id()])
          haveCommonPose = true;
      }
    }
    return haveCommonPose;
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::checkGraph()
  {
    double chi2 = 0.0;
    for (g2o::OptimizableGraph::EdgeSet::iterator it=_optimizer.edges().begin(); it!=_optimizer.edges().end(); it++){
      SLAMEdgeType* e1=dynamic_cast<SLAMEdgeType*>(*it);
      if(e1){
        chi2 = e1->chi2();
        if(chi2 < 0.0){
          std::cerr << "Error: NEGATIVE INFORMATION MATRIX... THIS SHOULD NOT HAPPEN..." << std::endl;
          std::cerr << "SLAMEdgeType " << chi2 << "\t";
          std::cerr << e1->information() << std::endl;
          e1->information().setIdentity();
          //exit(-1);
        }
        //cerr << e1->vertices()[0]->id() <<  " " << e1->vertices()[1]->id() << std::endl;
      } else {
        SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
        if(e2){
          chi2 = e2->chi2();
          if(chi2 < 0.0){
            std::cerr << "Error: NEGATIVE INFORMATION MATRIX... THIS SHOULD NOT HAPPEN..." << std::endl;
            std::cerr << "SensorEdgeType " << chi2 << "\t";
            e2->information().setIdentity();
            //std::cerr << e2->information() << std::endl;
            //exit(-1);
          }
//           if(e2->level() != 0)
//             std::cerr << PVAR(e2->level()) << std::endl;
          //cerr << e2->vertices()[0]->id() <<  " " << e2->vertices()[1]->id() << std::endl;
        } else {
          DataAssociationEdgeType* e3=dynamic_cast<DataAssociationEdgeType*>(*it);
          if(e3){
            chi2 = e3->chi2();
            if(chi2 < 0.0){
              std::cerr << "Error: NEGATIVE INFORMATION MATRIX... THIS SHOULD NOT HAPPEN..." << std::endl;
              std::cerr << "DataAssociationEdge " << chi2 << "\t";
              //std::cerr << e3->information() << std::endl;
              e3->information().setIdentity();
              //exit(-1);
            }
            //cerr << e3->vertices()[0]->id() <<  " " << e3->vertices()[1]->id() << std::endl;
          } else {
            std::cerr << "Warning: all edge casts went wrong..." << std::endl;
          }
        }
      }
    }

    for (g2o::OptimizableGraph::VertexIDMap::const_iterator it=_optimizer.vertices().begin(); it!=_optimizer.vertices().end(); it++){
      PoseVertex* v=dynamic_cast<PoseVertex*>(it->second);
      if(v){
        if(v->edges().size() == 0){
          std::cerr << "WARNING: ";
          std::cerr << PVAR(v->edges().size()) << std::endl;
        }
      } else {
        PointVertex* vp=dynamic_cast<PointVertex*>(it->second);
        if(vp)
          if(vp->edges().size() == 0){
            std::cerr << PVAR(vp->edges().size()) << std::endl;
          }
      }
    }
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  Observation<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex>& 
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::observation(int i)
  {
    return _verticies_observations[i];
  }

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  std::vector<typename SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::PointVertex* >&  
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::getPointVertices(int poseId, int level = 0)
  {

    if(_hierarchical_point_vertex_index.size() == 0){
      #pragma omp single
      for(size_t i = 0; i < _edges_observations.size(); ++i)
      {
        SensorEdgeType*& edge = _edges_observations[i];
        PointVertex* v=static_cast<PointVertex*>(edge->vertices()[1]);
        max_level_ = max(max_level_, edge->level());
        //std::cerr << i << " " << PVAR(edge->vertices()[0]->id()) << " " << PVAR(edge->level()) << std::endl;
        for(int j = edge->level(); j >= 0; --j)
          _hierarchical_point_vertex_index[edge->vertices()[0]->id()][j].push_back(v);
      }
    }
    return _hierarchical_point_vertex_index[poseId][level];
  }
    

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  void
  SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>::checkConnectivity()
  {
//     for(size_t i = 0; i < _verticies_points.size(); ++i){
//       PointVertex*& v = _verticies_points[i];
//       for (g2o::OptimizableGraph::EdgeSet::iterator it=v->edges().begin(); it!=v->edges().end(); it++){
//         SensorEdgeType* e2=dynamic_cast<SensorEdgeType*>(*it);
//         if(e2){
// //           if(!e2->vertices()[0]){
// //             std::cerr << __PRETTY_FUNCTION__ << " vertex not available. This SHOULD not happen! " << std::endl;
// //             exit(0);
// //           }
// //           parentMatch[e2->vertices()[0]->id()] = true;
//         }
//         DataAssociationEdgeType* e3=dynamic_cast<DataAssociationEdgeType*>(*it);
//         if(e3){
//           
//         }
//       }
//     }
  }

} //end namespace