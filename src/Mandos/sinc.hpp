#ifndef MANDOS_SINC_H_
#define MANDOS_SINC_H_

#include <Mandos/linear_algebra.hpp>

namespace mandos
{

Scalar sinc(Scalar x);

Vec3 grad_sinc(const Vec3& theta);

Mat3 hess_sinc(const Vec3& theta);

}  // namespace mandos

#endif  // MANDOS_SINC_H_
