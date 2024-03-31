#pragma once

#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>

namespace Flowtastic
{

using Vector2 = xt::xtensor_fixed<double, xt::xshape<2>>;
using Vector3 = xt::xtensor_fixed<double, xt::xshape<3>>;
using MatrixX = xt::xtensor<double, 2>;
using VectorX = xt::xtensor<double, 1>;

} // namespace Flowtastic