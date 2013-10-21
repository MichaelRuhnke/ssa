#ifndef STANDARD_CAMERA_H
#define STANDARD_CAMERA_H

namespace Vis {

class StandardCamera : public qglviewer::Camera
{
  public:
    StandardCamera() : qglviewer::Camera(), _standard(true), _closeRange(0.001), _farRange(1000.) {};

    float zNear() const {
      if (_standard) 
        return _closeRange; 
      else 
        return Camera::zNear(); 
    }

    float zFar() const
    {  
      if (_standard) 
        return _farRange; 
      else 
        return Camera::zFar();
    }

    void toggleMode() {_standard = !_standard;}
  
    void setFarRange(float r) { _farRange = r;}
    void setCloseRange(float r) { _closeRange = r;}

    bool isStandard() const {return _standard;}

  private:
    bool _standard;
    float _closeRange;
    float _farRange;
};

} // end namespace

#endif
