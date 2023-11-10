#ifndef GRAVITY_H_
#define GRAVITY_H_

struct GravityParamters {

};

struct Gravity {
  Gravity(unsigned int index, unsigned int dof, GravityParamters params)
    : index(index), dof(dof),paramters(params) {}

  const GravityParamters paramters;
  const unsigned int index;
  const unsigned int dof;
};

#endif // GRAVITY_H_
