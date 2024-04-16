#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include <cmath>

namespace Flowy::Math
{
inline constexpr double pi = 3.14159265358979323846264338327950288419716939937510;

inline double beta_pdf( double x, double alpha, double beta )
{
    double normalization = std::tgamma( alpha ) * std::tgamma( beta ) / std::tgamma( alpha + beta );
    return std::pow( x, alpha - 1.0 ) * std::pow( 1.0 - x, beta - 1.0 ) / normalization;
}

} // namespace Flowy::Math
