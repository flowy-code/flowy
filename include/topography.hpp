#pragma once
#include "asc_file.hpp"
#include "definitions.hpp"
#include <utility>

namespace Flowtastic
{

class Topography
{
public:
    Topography( const AscFile & asc_file )
            : height_data( asc_file.height_data ), x_data( asc_file.height_data ), y_data( asc_file.y_data ){};

    Topography() = default;

    MatrixX height_data{};
    VectorX x_data{};
    VectorX y_data{};

    inline double cell_size()
    {
        return x_data[1] - x_data[0];
    };

    // Calculate the height and the slope at coordinates
    // via linear interpolation from the square grid
    std::pair<double, Vector2> height_and_slope( const Vector2 & coordinates );
};

} // namespace Flowtastic