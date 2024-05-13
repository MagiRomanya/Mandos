#ifndef SINC_H_
#define SINC_H_

#include "linear_algebra.hpp"

Scalar sinc(Scalar x);

Vec3 grad_sinc(const Vec3& theta);

Mat3 hess_sinc(const Vec3& theta);

#endif // SINC_H_
