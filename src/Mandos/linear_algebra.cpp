#include <Mandos/linear_algebra.hpp>

namespace mandos
{

Eigen::Matrix<Scalar, 3, 9> vectorized_levi_civita()
{
    Eigen::Matrix<Scalar, 3, 9> e{{0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0},  //
                                  {0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},  //
                                  {0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0}};
    return e;
}
}  // namespace mandos