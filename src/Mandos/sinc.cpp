#include <Mandos/sinc.hpp>

namespace
{
constexpr mandos::Scalar SINC_THRESHOLD = 1e-6;
}

namespace mandos
{

Scalar sinc(Scalar x)
{
    if (abs(x) < SINC_THRESHOLD) {
        const Scalar x2 = x * x;
        const Scalar x4 = x2 * x2;
        return 1.0 - x2 / 6.0 + x4 / 120.0;
    }
    return std::sin(x) / x;
}

inline Scalar sinc_g_function(Scalar x)
{
    if (abs(x) < SINC_THRESHOLD) {
        const Scalar x2 = x * x;
        const Scalar x4 = x2 * x2;
        return -1.0 / 3.0 + x2 / 30.0 + x4 / 840.0;
    }

    return (x * std::cos(x) - std::sin(x)) / (x * x * x);
}

inline Scalar sinc_h_function(Scalar x)
{
    const Scalar x2 = x * x;
    const Scalar x4 = x2 * x2;
    if (abs(x) < SINC_THRESHOLD) {
        return 1.0 / 15.0 + x2 / 210.0 + x4 / 7560.0;
    }

    return (-x2 * std::sin(x) - 3.0 * x * std::cos(x) + 3.0 * std::sin(x)) / (x4 * x);
}

Vec3 grad_sinc(const Vec3& theta)
{
    const Scalar x = theta.norm();
    const Scalar g = sinc_g_function(x);
    return g * theta;
}

Mat3 hess_sinc(const Vec3& theta)
{
    const Scalar x = theta.norm();
    const Scalar g = sinc_g_function(x);
    const Scalar h = sinc_h_function(x);
    return h * theta * theta.transpose() + g * Mat3::Identity();
}

}  // namespace mandos
