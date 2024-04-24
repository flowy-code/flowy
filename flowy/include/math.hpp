#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include <fmt/format.h>
#include <cmath>
#include <stdexcept>
#include <string>

namespace Flowy::Math
{
inline constexpr double pi = 3.14159265358979323846264338327950288419716939937510;

inline double beta_pdf( double x, double alpha, double beta )
{
    double normalization = std::tgamma( alpha ) * std::tgamma( beta ) / std::tgamma( alpha + beta );
    return std::pow( x, alpha - 1.0 ) * std::pow( 1.0 - x, beta - 1.0 ) / normalization;
}

#define FLOWY_CHECK( x ) Flowy::Math::check_number( #x, x, __FILE__, __LINE__ )

template<typename N>
inline void check_number( const std::string & name, const N & var, const std::string & file, int line )
{
    if( !std::isfinite( var ) || std::isnan( var ) || std::isinf( var ) )
    {
        throw std::runtime_error(
            fmt::format( "File {}:{} - '{}' is {}. This should NEVER happen.", file, line, name, var ) );
    }
}

} // namespace Flowy::Math
