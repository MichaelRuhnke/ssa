// Sparse Surface Optimization
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

#include <iostream>
#include <signal.h>
#include <limits>

#include "g2o/stuff/macros.h"
#include "ssa/core/ssa_graph_3d.h"
#include "ssa/data_association/data_association.h"


using namespace std;
using namespace ssa;

void evaluateResult(NearestNeighborBruteForceT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::CorrespondenceList& correspondences, std::map<int, std::map<int, std::map<int,  NearestNeighborBruteForceT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::Correspondence* > > >& ground_truth){
  int possibleMatches = 0;
  int matches = 0;
  int distMatches = 0;
  int groundTruthDistMatches = 0;
  int kDTreeBetter = 0;
  int kDTreeWorse = 0;
  for(size_t i = 0; i < correspondences.size(); ++i)
  {
    NearestNeighborBruteForceT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::Correspondence* correspondence = correspondences[i];
    NearestNeighborBruteForceT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::Correspondence* reference_correspondence =     ground_truth[correspondence->query_point->parentVertex()->id()][correspondence->cor_point->parentVertex()->id()][correspondence->query_point->id()];

    if(correspondence && ! reference_correspondence)
    {
       std::cerr << "Warning no reference found for correspondence "  << (correspondence->query_point->id()) << " -> " << (correspondence->cor_point->id()) << " " << (correspondence->query_point->estimate() - correspondence->cor_point->estimate()).squaredNorm() << std::endl;
    }

    if(correspondence && reference_correspondence){
      possibleMatches++;
      bool almostSimilarDistance = fabs(correspondence->sqrDistance - reference_correspondence->sqrDistance) < std::numeric_limits<float>::epsilon();
      ///check if kdTree has better or worse solution
      if(correspondence->sqrDistance < reference_correspondence->sqrDistance && !almostSimilarDistance)
      {
        kDTreeBetter++;
      } else {
        if(correspondence->sqrDistance > reference_correspondence->sqrDistance && !almostSimilarDistance)
          kDTreeWorse++;
//         if(almostSimilarDistance && correspondence->cor_point->id() != reference_correspondence->cor_point->id())
//           std::cerr << PVAR(correspondence->cor_point->id()) << " " << PVAR(reference_correspondence->cor_point->id()) << std::endl;
      }
        
      if(correspondence->cor_point->id() == reference_correspondence->cor_point->id() && correspondence->query_point->id() == reference_correspondence->query_point->id()){
        matches++;

        float trueDistance = (correspondence->cor_point->estimate() - correspondence->query_point->estimate()).squaredNorm();
        if (fabs(correspondence->sqrDistance - trueDistance) < std::numeric_limits<float>::epsilon())
          distMatches++;
        if (fabs(reference_correspondence->sqrDistance - trueDistance) < std::numeric_limits<float>::epsilon())
          groundTruthDistMatches++;
      } else {
//         std::cerr << "ERROR: " << (correspondence->query_point->id()) << " -> " <<(correspondence->cor_point->id()) << "(" << correspondence->sqrDistance << ")" << " vs. reference " << (reference_correspondence->query_point->id()) << " -> " << (reference_correspondence->cor_point->id())  << "(" << reference_correspondence->sqrDistance << ")" << std::endl;
      }
    }
  }
  cerr << "Evaluation: " << std::endl;
  cerr << "Match percentage   \t " << ((float) possibleMatches / (float) correspondences.size()) * 100.0f <<"% " << std::endl;  
  cerr << "Id matches         \t " << ((float) matches / (float) possibleMatches) * 100.0f <<"% " << std::endl;
  cerr << "result better dist \t " << ((float) kDTreeBetter / (float) possibleMatches) * 100.0f <<"% " << std::endl;
  cerr << "result worse dist  \t " << ((float) kDTreeWorse / (float) possibleMatches) * 100.0f <<"% " << std::endl;
  cerr << "Distance matches   \t " << ((float) distMatches / (float) possibleMatches) * 100.0f <<"% " << std::endl;
  cerr << "GT distance matches\t " << ((float) groundTruthDistMatches / (float) possibleMatches) * 100.0f <<"% " << std::endl;
}


const char *message[]={
  "ssa_test_data_association: tests and compares all available data association strategies",
  "usage ssa_test_data_association [options] <ssa_file>",
  "options:",
  "-ini    provide configuration.",
  0
};

int main(int argc, char **argv)
{

  const char**v=message;
  while (*v){
    cout << *v << endl;
    v++;
  }

  // HACK reset numeric locale, we need C
  setlocale (LC_NUMERIC,"C");
  int c=1;
  const char* logfile=0;
  bool useConfigFile = false;
  const char* configFile=0;

  while (c<argc){
    if (!strcmp(argv[c],"-debug")){
//       debug=true;
      c++;
    } else 
    if (!strcmp(argv[c],"-ini")){
      useConfigFile=true;
      c++;
      configFile = argv[c];
      c++;
    } else 
    if (! logfile){
      logfile=argv[c];
      c++;
      break;
    }
  }

  if(!logfile){
    cerr << "Error: no input ssa file given." << endl;
    exit(0);
  }

  double timing = get_time();
  /** loading graph  */
  std::cerr << "loading ssa graph                 \t ";
  SparseSurfaceAdjustmentGraph3D ssaGraph;
  if(logfile){
    ssaGraph.load(logfile);
  }
  std::cerr << " [done] " << (get_time() - timing) * 1000.0 << "ms." << std::endl;

  timing = get_time();
  std::cerr << "loading params from file          \t ";
  SparseSurfaceAdjustmentParams params;
  /** Load configuration from ini file */
  if(useConfigFile){
    ifstream configStream(configFile);
    params.readParams(configStream);
    configStream.close();
    //params.printParams();
  }
  std::cerr << " [done] " << (get_time() - timing) * 1000.0 << "ms." << std::endl;

  timing = get_time();
  std::cerr << "creating neighbor cache           \t ";
  ssaGraph.fillNeighborCache(params);
  std::cerr << " [done] " << (get_time() - timing) * 1000.0 << "ms." << std::endl; 

  timing = get_time();
  std::cerr << "Brute force (overlap heuristic)   \t ";
  DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> da_brute_force;
  da_brute_force.setStrategy(DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::BRUTEFORCE);
  da_brute_force.applyBasedOnOverlap(ssaGraph, params, 0);
  std::cerr << " [done] " << da_brute_force.correspondences.size() << " results in " << (get_time() - timing) * 1000.0 << "ms." << std::endl;
  da_brute_force.correspondences.clear();

  timing = get_time();
  std::cerr << "Brute force correspondence search \t ";
  da_brute_force.apply(ssaGraph, params, 0);
  std::cerr << " [done] " << da_brute_force.correspondences.size() << " results in " << (get_time() - timing) * 1000.0 << "ms." << std::endl;

  std::map<int, std::map<int, std::map<int,  NearestNeighborBruteForceT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::Correspondence* > > > ground_truth;
  for(size_t i = 0; i < da_brute_force.correspondences.size(); ++i)
  {
    NearestNeighborBruteForceT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::Correspondence* correspondence = da_brute_force.correspondences[i];
    ground_truth[correspondence->query_point->parentVertex()->id()][correspondence->cor_point->parentVertex()->id()][correspondence->query_point->id()] = correspondence;
  }

  timing = get_time();
  std::cerr << "kDTree correspondence search      \t ";
  DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> da_kdtree;
  da_kdtree.setStrategy(DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::KDTREE);
  da_kdtree.apply(ssaGraph, params, 0);
  std::cerr << " [done] " << da_kdtree.correspondences.size() << " results in " << (get_time() - timing) * 1000.0 << "ms." << std::endl;
  evaluateResult(da_kdtree.correspondences, ground_truth);

  /** retest with overlap heuristic */
  da_kdtree.correspondences.clear();
  timing = get_time();
  std::cerr << "kDTree (overlap heuristic)       \t ";
  da_kdtree.applyBasedOnOverlap(ssaGraph, params, 0);
  std::cerr << " [done] " << da_kdtree.correspondences.size() << " results in " << (get_time() - timing) * 1000.0 << "ms." << std::endl;
  evaluateResult(da_kdtree.correspondences, ground_truth);


  timing = get_time();
  std::cerr << "local with sparse kDTree search   \t ";
  DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> da_sparse_kdtree_local;
  da_sparse_kdtree_local.setStrategy(DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::SPARSE_KDTREE_LOCAL);
  da_sparse_kdtree_local.apply(ssaGraph, params, 0);
  std::cerr << " [done] " << da_sparse_kdtree_local.correspondences.size() << " results in " << (get_time() - timing) * 1000.0 << "ms." << std::endl;

  evaluateResult(da_sparse_kdtree_local.correspondences, ground_truth);


//   timing = get_time();
//   cerr << "calculation overlap matrix on highest level: \t" << std::endl;
//   da_kdtree.overlap(ssaGraph, params, 0);
//   std::cerr << " [done] " << " in " << (get_time() - timing) * 1000.0 << "ms." << std::endl;
}
