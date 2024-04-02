#pragma once
#include "asc_file.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include <utility>

namespace Flowtastic
{

class Topography
{
public:
    // The indices are all inclusive, i.e.
    // x limits = [idx_x_lower, idx_x_higher]
    // y limits = [idx_y_lower, idx_y_higher]
    struct BoundingBox
    {
        int idx_x_lower{};
        int idx_x_higher{};
        int idx_y_lower{};
        int idx_y_higher{};
    };

    Topography( const AscFile & asc_file )
            : height_data( asc_file.height_data ), x_data( asc_file.x_data ), y_data( asc_file.y_data ){};

    Topography( const MatrixX & height_data, const VectorX & x_data, const VectorX & y_data )
            : height_data( height_data ), x_data( x_data ), y_data( y_data ){};

    Topography() = default;

    MatrixX height_data{}; // The heights of the cells
    VectorX x_data{};
    VectorX y_data{};

    inline double cell_size()
    {
        return x_data[1] - x_data[0];
    };

    // Calculate the height and the slope at coordinates
    // via linear interpolation from the square grid
    std::pair<double, Vector2> height_and_slope( const Vector2 & coordinates );

    // Compute the indices of a rectangular bounding box
    // The box is computed such that a circle with centered at 'center' with radius 'radius'
    // Is completely contained in the bounding box
    BoundingBox bounding_box( const Vector2 & center, double radius );

    // Find all the cells that intersect the lobe
    std::vector<std::array<int, 2>> get_cells_intersecting_lobe( const Lobe & lobe );

    // Find the fraction of the cells covered by the lobe by rasterizing each cell
    // into a grid of N*N points
    // This returns a vector of pairs
    // - the first entry of each pair contains an array<int, 2> with the idx_i, idx_j of the intersected cell
    // - the second entry contains the fraction of the cell that is covered by the ellips
    std::vector<std::pair<std::array<int, 2>, double>> compute_intersection( const Lobe & lobe, int N = 15 );

    // Adds the lobe thickness to the topography, according to its fractional intersection with the cells
    void add_lobe( const Lobe & lobe );

    // Figure out which cell a given point is in, returning the indices of the lowest left corner
    std::array<int, 2> locate_point( const Vector2 & coordinates );
};

} // namespace Flowtastic