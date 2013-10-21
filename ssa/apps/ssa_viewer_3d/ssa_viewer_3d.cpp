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

#include <iostream>
#include <signal.h>
#include <QApplication>
#include <QtGui/QMainWindow>
#include "ui_base_main_window.h"

#include "ssa/core/allocate_solver.h"
#include "ssa/core/ssa_graph_3d.h"
#include "ssa/core/sparse_surface_adjustment.h"
#include "ssa/data_association/local_refinement.h"
#include "ssa/data_association/representative_subset.h"
//#include "ssa/solvers/pcg_cuda/linear_solver_pcg_cuda.h"
// #include "ssa/solvers/pcg_cusp/linear_solver_pcg_cusp.h"

using namespace std;
using namespace ssa;

static bool running = true;
unsigned int breakCounter = 0;

void sighandler(int sig)
{
    cout << endl << "Signal " << sig << " caught..." << endl;

   running = false;
   breakCounter++;
   if(breakCounter >= 3){
      exit(1);
   }
}


int main(int argc, char **argv)
{
  signal(SIGINT, &sighandler);
  QApplication qapp(argc, argv);
  // HACK reset numeric locale, we need C
  setlocale (LC_NUMERIC,"C");
  //bool debug = false;
  bool restoreViewerState = false;
  bool dumpScreenshots = false;
  bool saveOutput = false;
  char *outfile = 0;

  const char* logfile=0;
  const char* viewerStateFile=0;
  bool useConfigFile = false;
  const char* configFile=0;
  int level = 1;
  int c=1;
  while (c<argc){
    if (!strcmp(argv[c],"-cam")){
      restoreViewerState=true;
      c++;
      viewerStateFile = argv[c];
      c++;
    } else
    if (!strcmp(argv[c],"-level")){
      c++;
      level = atoi(argv[c]);
      c++;
    } else
    if (!strcmp(argv[c],"-ini")){
      useConfigFile=true;
      c++;
      configFile = argv[c];
      c++;
    } else
    if (!strcmp(argv[c],"-video")){
      dumpScreenshots=true;
      c++;
      viewerStateFile = argv[c];
      c++;
    } else
    if (!strcmp(argv[c],"-save")){
      saveOutput=true;
      c++;
      outfile = argv[c];
      c++;
    } else
    if (! logfile){
      logfile=argv[c];
      c++;
      break;
    }
  }

  //initialize optimizer
  SparseSurfaceAdjustmentT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> ssa;
  ssa.setVerbose(true);
  g2o::BlockSolverX::LinearSolverType* linearSolver = AllocateLinearSolver<g2o::BlockSolverX>(1); //CHOLMOD
  ssa.setSolver(linearSolver);

  if(useConfigFile){
    ifstream configStream(configFile);
    ssa.params().readParams(configStream);
    configStream.close();
    ssa.params().printParams();
  }

  if (! logfile){

  } else {
    ssa.graph()->load(logfile);
  }

  ssa.graph()->fillNeighborCache(ssa.params());

    QMainWindow* mw = new QMainWindow;
    Ui_MainWindow umw;
    umw.setupUi(mw);

    umw.ssaGLWidget->setMap(ssa.graph());
    umw.spinBox_2->setValue(level);
    mw->show();
    umw.ssaGLWidget->setBackgroundColor(qRgb(211, 211, 211));
    if(restoreViewerState){
      umw.ssaGLWidget->setStateFileName(viewerStateFile);
      umw.ssaGLWidget->restoreStateFromFile();
    }

    umw.ssaGLWidget->_iterations = ssa.params().ssaIterations;
    umw.doubleSpinBox_3->setValue(ssa.params().normalShooting.stepSize);
    umw.doubleSpinBox_4->setValue(ssa.params().normalShooting.steps);
    umw.doubleSpinBox_5->setValue(ssa.params().normalShooting.maxAngleDifference);

    QObject::connect( &(ssa.graph()->_optimizer), SIGNAL(iterationDone()), umw.ssaGLWidget, SLOT(Gen3DObjectList_update()));
    QObject::connect( &(ssa.graph()->_optimizer), SIGNAL(iterationDone()), umw.ssaGLWidget, SLOT(updateGL()));

    if(dumpScreenshots)
      QObject::connect( &(ssa.graph()->_optimizer), SIGNAL(iterationDone()), umw.ssaGLWidget, SLOT(saveSnapshotVideo()));
    bool needReDraw = true;

    if(ssa.graph()->_verticies_poses.size() > 0){
      g2o::VertexSE3* v = ssa.graph()->_verticies_poses[0];

      qglviewer::Vec initialNodePose = qglviewer::Vec(v->estimate().translation()(0), v->estimate().translation()(1), v->estimate().translation()(2)+20.0);
      umw.ssaGLWidget->camera()->setPosition(initialNodePose);

      initialNodePose = qglviewer::Vec(v->estimate().translation()(0), v->estimate().translation()(1), v->estimate().translation()(2));
      umw.ssaGLWidget->camera()->setSceneCenter(initialNodePose);

      umw.ssaGLWidget->camera()->setUpVector(qglviewer::Vec(0.0, 0.0, 1.0));
      umw.ssaGLWidget->camera()->lookAt(initialNodePose);
    }
    cerr << ssa.graph()->_verticies_poses.size() << endl;
    ///rebuild / build display lists
    umw.ssaGLWidget->Gen3DObjectList_update();

    while(mw->isVisible() && running) {
       qapp.processEvents();
       if(umw.ssaGLWidget->_optimize){
          umw.ssaGLWidget->_optimize = false;

        if(dumpScreenshots)
          umw.ssaGLWidget->saveSnapshotVideo();

        umw.stateLabel->setText(QString("running"));
        QPalette palette;
        palette.setColor(umw.stateLabel->backgroundRole(), Qt::yellow);
        umw.stateLabel->setPalette(palette);
        umw.stateLabel->setAutoFillBackground(true);
        mw->repaint();
        qapp.processEvents();

        //filling selected gui params into ssa params
        ssa.params().ssaIterations = umw.ssaGLWidget->_iterations;
        //ssa.params().g2oIterations = 6;
        ssa.params().normalShooting.stepSize = umw.doubleSpinBox_3->value();
        ssa.params().normalShooting.steps = (int) umw.doubleSpinBox_4->value();
        ssa.params().normalShooting.maxAngleDifference = umw.doubleSpinBox_5->value();
        ssa.params().printParams();
        cerr << "optimizing level " << umw.spinBox_2->value() << endl;
        ssa.optimize(umw.spinBox_2->value());
        //ssa.optimize(ssa.graph()->getMaxLevel());
        umw.ssaGLWidget->Gen3DObjectList_update();
        umw.ssaGLWidget->updateGL();

        umw.stateLabel->setText(QString("ready"));
        palette.setColor(umw.stateLabel->backgroundRole(), Qt::green);
        umw.stateLabel->setPalette(palette);

        umw.stateLabel->setAutoFillBackground(true);
        needReDraw = true;

        if(dumpScreenshots){
          umw.ssaGLWidget->Gen3DObjectList_update();
          umw.ssaGLWidget->updateGL();
          needReDraw = false;
          umw.ssaGLWidget->saveSnapshotVideo();
        }

        if(saveOutput){
          cerr << "wrote resulting graph to " << outfile << endl;
          ssa.graph()->save(outfile);
        }


//         cerr << "calculation took " << get_time() - startTime << " seconds." << endl;
       }

       if(umw.ssaGLWidget->_solveDA){
          umw.ssaGLWidget->_solveDA = false;

	  std::cerr << "create subset \t";
	  RepresentativeSubsetT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> subset;
	  subset.createSubset(ssa.graph(), ssa.params());
	  umw.ssaGLWidget->Gen3DObjectList_update();
          umw.ssaGLWidget->updateGL();
	  std::cerr << "done" << std::endl;;

//         ssa.graph()->calcMeanCov(ssa.params());
//         ssa.graph()->dropDataAssociation();
//         DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov> da_kdtree;
//         da_kdtree.setStrategy(DataAssociationT<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::KDTREE);
//         da_kdtree.apply(*ssa.graph(), ssa.params(), 0);
//
//         std::cerr << "Filter outlier...";
//         ssa.graph()->filterOutlier(ssa.params().outlierRejectionMinConnectedNeighbors);
//         std::cerr << "done" << std::endl;
//         std::cerr << "optimize color...";
//         ssa.optimizeColorsIncidenteWeight();
//         std::cerr << "done" << std::endl;
//         umw.ssaGLWidget->Gen3DObjectList_update();
//         ssa.graph()->dropDataAssociation();

//           if(ssa.graph()->_edges_data_association.size() == 0){
//
//             //filling selected gui params into ssa params
//             ssa.params().ssaIterations = umw.ssaGLWidget->_iterations;
//             ssa.params().normalShooting.stepSize = umw.doubleSpinBox_3->value();
//             ssa.params().normalShooting.steps = (int) umw.doubleSpinBox_4->value();
//             ssa.params().normalShooting.maxAngleDifference = umw.doubleSpinBox_5->value();
//             ssa.params().printParams();
//
//             ssa.graph()->dropDataAssociation();
//             //ssa.graph()->calcMeanCov(ssa.params().normalExtractionMaxNeighborDistance);
//             //NormalShootingFlann<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::shootNormals(*ssa.graph(), ssa.params().normalShooting);
//             double timing = get_time();
//   //           NearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::apply(*ssa.graph(), ssa.params().nearestNeighbor, umw.spinBox_2->value());
//             NearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::apply(*ssa.graph(), ssa.params().nearestNeighbor, umw.spinBox_2->value());
//             std::cerr << "created " <<  ssa.graph()->_edges_data_association.size() << " d.a. edges in " << (get_time() - timing) * 1000 << "ms."<< std::endl;
//           } else {
//             double timing = get_time();
//             LocalRefinement<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::apply(*ssa.graph(), ssa.params(), umw.spinBox_2->value());
//             std::cerr << "updated " <<  ssa.graph()->_edges_data_association.size() << " d.a. edges in " << (get_time() - timing) * 1000 << "ms."<< std::endl;
//           }

//           ssa.graph()->dropDataAssociation();
//           SparseNearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::apply(*ssa.graph(), ssa.params(), umw.spinBox_2->value());
          needReDraw = true;
       }

       if(umw.ssaGLWidget->_prune){
          umw.ssaGLWidget->_prune = false;
          double timing = ssa::get_time();
          //filling selected gui params into ssa params
          //NearestNeighbor<g2o::EdgeSE3, ssa::EdgeSE3PointXYZCov, ssa::EdgePointXYZCovPointXYZCov>::pruneToMinSpanTree(*ssa.graph());
          ssa.graph()->pruneGraph(ssa.params().normalShooting.stepSize);
          std::cerr << "pruneGraph took " << ssa::get_time() - timing << " s." << std::endl;
          needReDraw = true;
       }

       if(needReDraw){
         umw.ssaGLWidget->updateGL();
         needReDraw = false;
       }
       usleep(1000);
    }
  return 0;
}
