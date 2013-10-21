#include "ssa/core/ssa_graph_3d.h"
#include <GL/gl.h>
#include <GL/glut.h>
#include <QGLViewer/qglviewer.h>

namespace ssa {
  using namespace Eigen;
  using namespace std;

class SSARGBAglWidget: public QGLViewer
{
  Q_OBJECT
  public:
    SSARGBAglWidget();
    SSARGBAglWidget(QWidget*& widget);
    ~SSARGBAglWidget();

  inline void setInitialParams(){
    _backgroundHeight = -0.01;
    _height = 0.01;
    _landmarkSize = 0.01;
    _poseSize = 0.2;
    _optimize = false;
    _iterations = 10;
    _DAmaxdist = 0.2;
    _pointSize = 0.01;
    _drawVertices = true;
    _drawLandmarkVertices = true;
    _drawNormals = true;
    _drawCorrespondences = true;
    _drawMesh = true;
    _drawGrid = false;
    _prune = false;
    _solveDA = false;
    setShortcut(EXIT_VIEWER, Qt::Key_Q);
    setSnapshotCounter(0);
    _level = 0;
  }

  void init();

  void draw();
  virtual void drawWithNames();
  void drawEdges();
  void drawCovariance();
  void drawMeshEdges();
  void drawVertices();

  /** draws circle with radius r at the current pose
  */
  void drawEllipsoid(double& l1, double& l2, double& l3);
  void drawCircle(GLfloat radius, int segments);

  void drawPoseCircle(double& x, double& y, double& theta);
  void drawPoseCircleWithName(double& x, double& y, double& z);

  void drawNode(Vector2d& pose);
  void drawGrid();
  void drawPoseBox();
  void drawBox(GLfloat l, GLfloat w, GLfloat h);

  void postSelection(const QPoint& point);

  inline void setMap(SparseSurfaceAdjustmentGraphRGBA* accurateMap){
    _ssa_graph = accurateMap;
    Gen3DObjectList_update();
  }

  GLint Gen3DObjectList_poseVertices();
  GLint Gen3DObjectList_obsVertices();
  GLint Gen3DObjectList_normals();
  GLint Gen3DObjectList_mesh();

  public slots:
      void setOptimize();
      inline void prune() {_prune = true; };
      inline void solveDA() {_solveDA = true; };
      void openFile();
      void saveFile();
      inline void setIterations(int it){ _iterations = it;};
      inline void setDAMaxDist(double i){ _DAmaxdist = i;};
      inline void setPointSize(double i){ _pointSize = i; updateGL(); cerr << "pointSize " << _pointSize << endl; Gen3DObjectList_update();};
      inline void setDrawGrid(int i){ _drawGrid = (bool) i; updateGL();}
      inline void setDrawVertices(int i){ _drawVertices = (bool) i; updateGL();}
      inline void setDrawLandmarkVertices(int i){ _drawLandmarkVertices = (bool) i; updateGL();}
      inline void setDrawNormals(int i){ _drawNormals = (bool) i; updateGL();}
      inline void setDrawCorrespondences(int i){ _drawCorrespondences = (bool) i; updateGL();}
      inline void setDrawMesh(int i){ _drawMesh = (bool) i; updateGL();}
      inline void setStrategie(int i){ std::cerr << __PRETTY_FUNCTION__ << " not implemented!" << i << std::endl;}
      inline void saveFMap(){ std::cerr << __PRETTY_FUNCTION__ << " not implemented!" << std::endl;}
      inline void setDrawLevel(int i){ _level = i;};
      void saveSnapshotCustom();
      void saveSnapshotVideo();
      void Gen3DObjectList_update();
      void Gen3DObjectList_updateSelection();

  public:
  QWidget*               _parent;
  SparseSurfaceAdjustmentGraphRGBA*  _ssa_graph;

  bool _optimize;
  bool _drawGrid;
  bool _drawVertices;
  bool _drawLandmarkVertices;
  bool _drawNormals;
  bool _drawCorrespondences;
  bool _drawMesh;
  bool _prune;
  bool _solveDA;
  GLint _poseVertices;
  GLint _obsVertices;
  GLint _normals;
  GLint _mesh;

  unsigned int _iterations;
  double _DAmaxdist;

  double _backgroundHeight;
  double _height;
  double _landmarkSize;
  double _poseSize;
  double _pointSize;
  qglviewer::Vec _selectedPoint;
  int     _level;
};

}
