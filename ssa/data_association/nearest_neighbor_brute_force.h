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

#ifndef __SSA_DA_BRUTE_FORCE_T__
#define __SSA_DA_BRUTE_FORCE_T__

#include <vector>
#include <tr1/unordered_map> 
#include "ssa/core/ssa_graph.h"
#include "ssa/core/parameter.h"
#include "correspondence.h"
#include "correspondence_rejection.h"
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() 0
    #define omp_get_thread_num() 0
#endif

namespace ssa {

  template <typename EdgeType1, typename EdgeType2, typename EdgeType3>
  class NearestNeighborBruteForceT{

    typedef EdgeType1 SLAMEdgeType;
    typedef EdgeType2 SensorEdgeType;
    typedef EdgeType3 DataAssociationEdgeType;

    typedef typename SensorEdgeType::VertexXjType   PointVertex;

    public:
    typedef typename std::map<int, std::map<int, CorrespondenceT<PointVertex, DataAssociationEdgeType> > > CorrespondenceIdMap;
    typedef CorrespondenceT<PointVertex, DataAssociationEdgeType>                                          Correspondence;
    typedef typename std::deque< Correspondence*  >                                                        CorrespondenceList;
    typedef std::tr1::unordered_map<int, std::tr1::unordered_map<int, float > >                            OverlapMap;
    typedef std::vector< std::pair <int,int> >                                                             ScanPairVector;

    NearestNeighborBruteForceT();
    ~NearestNeighborBruteForceT();

    /** \brief creates data association edges for all available scan pairs  */
    void apply(CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    /** \brief creates data association edges for all available scan pairs  */
    void apply(OverlapMap overlap, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    /** \brief creates data association edges for all available scan pairs  */
    virtual void apply(ScanPairVector scanPairs, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, SparseSurfaceAdjustmentParams& params, int level);

    /** \brief creates data association edges for all available scan pairs  */
    void applyForScanPair(int id1, int id2, CorrespondenceList& resultingCorrespondences, SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, double sqrDistance, int level);

    /** \brief returns a vector with all scan pairs (given N scans this method will return (N*N-1) scan pairs!!)*/
    ScanPairVector getScanPairs(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph);

    /** \brief returns a vector with all incremental scan pairs (given N scans this method will return (N-1) scan pairs!!)*/
    ScanPairVector getIncrementalScanPairs(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph);

    /** \brief returns a vector with all scan pairs (based on a given overlap treshold)*/
    ScanPairVector getScanPairs(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, OverlapMap overlap, float minOverlap);

    void assign(CorrespondenceList& resultingCorrespondences, PointVertex*& reference, PointVertex*& correspondence, float distance);

    /** \brief creates an edge between reference and correspondence */
    void assign(SparseSurfaceAdjustmentGraphT<EdgeType1, EdgeType2, EdgeType3>& graph, PointVertex*& reference, PointVertex*& correspondence);

   };

   #include "nearest_neighbor_brute_force.hpp"
}

#endif
