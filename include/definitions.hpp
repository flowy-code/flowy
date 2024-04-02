#pragma once

#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#include "xtensor-blas/xlinalg.hpp"
#pragma GCC diagnostic pop

namespace Flowtastic
{

using Vector2 = xt::xtensor_fixed<double, xt::xshape<2>>;
using Vector3 = xt::xtensor_fixed<double, xt::xshape<3>>;
using MatrixX = xt::xtensor<double, 2>;
using VectorX = xt::xtensor<double, 1>;

} // namespace Flowtastic