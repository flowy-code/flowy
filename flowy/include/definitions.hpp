#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers

#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmanipulation.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <xtensor/xcsv.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace Flowy
{

using Vector2 = xt::xtensor_fixed<double, xt::xshape<2>>;
using Vector3 = xt::xtensor_fixed<double, xt::xshape<3>>;
using MatrixX = xt::xtensor<double, 2>;
using VectorX = xt::xtensor<double, 1>;

// The default NO_DATA_VALUE to be used for elevation data
constexpr double DEFAULT_NO_DATA_VALUE_HEIGHT = -9999;

// The default NO_DATA_VALUE to be used for thickness data
constexpr double DEFAULT_NO_DATA_VALUE_THICKNESS = 0.0;

enum StorageDataType
{
    Short,
    Float,
    Double
};

} // namespace Flowy
