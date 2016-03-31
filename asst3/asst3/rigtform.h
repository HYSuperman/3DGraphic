#ifndef RIGTFORM_H
#define RIGTFORM_H

#include <iostream>
#include <cassert>

#include "matrix4.h"
#include "quat.h"

class RigTForm {
  Cvec3 t_; // translation component
  Quat r_;  // rotation component represented as a quaternion

public:
  RigTForm() : t_(0) {
    assert(norm2(Quat(1,0,0,0) - r_) < CS175_EPS2);
  }

  RigTForm(const Cvec3& t, const Quat& r): t_(t), r_(r){
    // asst3
  }

  explicit RigTForm(const Cvec3& t): t_(t), r_(Quat()){
    // asst3
  }

  explicit RigTForm(const Quat& r): t_(Cvec3()), r_(r) {
    // asst3
  }

  Cvec3 getTranslation() const {
    return t_;
  }

  Quat getRotation() const {
    return r_;
  }

  RigTForm& setTranslation(const Cvec3& t) {
    t_ = t;
    return *this;
  }

  RigTForm& setRotation(const Quat& r) {
    r_ = r;
    return *this;
  }

  Cvec4 operator * (const Cvec4& a) const {
    Cvec4 t = Cvec4(t_);
    Quat r = getRotation();
    return r*a + t;
    // asst3
  }

  RigTForm operator * (const RigTForm& a) const {
    Cvec4 t1 = Cvec4(getTranslation());
    Cvec4 t2 = Cvec4(a.getTranslation());
    Quat r1 = getRotation();
    Quat r2 = a.getRotation();
    return RigTForm(Cvec3(t1+r1*t2), r1*r2);
    // asst3
  }

};

inline RigTForm inv(const RigTForm& tform) {
  Quat r = tform.getRotation();
  Cvec4 t = Cvec4(tform.getTranslation());
  return RigTForm(-Cvec3(inv(r)*t), inv(r));
  // asst3
}

inline RigTForm transFact(const RigTForm& tform) {
  return RigTForm(tform.getTranslation());
}

inline RigTForm linFact(const RigTForm& tform) {
  return RigTForm(tform.getRotation());
}

inline Matrix4 rigTFormToMatrix(const RigTForm& tform) {
  Matrix4 T = Matrix4::makeTranslation(tform.getTranslation());
  Matrix4 R = quatToMatrix(tform.getRotation());
  return T * R;
  // asst3
}

// ---------------------------------------------
// Take the object frame and eye frame
// return (O)t(E)r
// ---------------------------------------------
inline RigTForm makeMixedFrame(const RigTForm& obj, const RigTForm& eye){
  // multiply transforlation of O and linearTrans of E.
  RigTForm r = transFact(obj)*linFact(eye);
  return r;

}

// --------------------------------------------
// do Q to O with respect to A
// --------------------------------------------
inline RigTForm doQtoOwrtA(const RigTForm& q, const RigTForm& o, const RigTForm& a){
  // just calculate...
  RigTForm r = a * q * inv(a) * o;
  return r;
}

#endif
