#pragma once
#include "asc_file.hpp"
#include "definitions.hpp"

namespace Flowtastic
{

class Topography
{
public:
    Topography( const AscFile & asc_file )
            : height_data( asc_file.height_data ), x_data( asc_file.height_data ), y_data( asc_file.y_data ), cell_size( asc_file.cell_size ){};
    
    Topography() = default;

    MatrixX height_data{};
    VectorX x_data{};
    VectorX y_data{};
    double cell_size{};

    // Calculate the slope at coordinates, via linear interpolation from the square grid
    Vector3 slope( const Vector2 & coordinates );
};

} // namespace Flowtastic