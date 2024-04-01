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
    /*
    The indices are all inclusive, i.e.
    x limits = [idx_x_lower, idx_x_higher]
    y limits = [idx_y_lower, idx_y_higher]
    */
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

    BoundingBox bounding_box( const Vector2 & center, double radius );

    // Test if a given point lies within the cell
    // Note that the left and bottom border of the cell are included, while the right and top are not
    bool point_in_cell( int idx_i, int idx_j, const Vector2 & point );

    // Test if a line given by y(x) = slope_xy * x + offset intersects the cell
    bool line_intersects_cell(
        int idx_i, int idx_j, double slope_xy, double offset, double x_min = -std::numeric_limits<double>::infinity(),
        double x_max = std::numeric_limits<double>::infinity() );

    // Test if a line segment between x1 and x2 intersects the cell
    bool line_segment_intersects_cell( int idx_x, int idx_y, Vector2 x1, Vector2 x2 );

    // Find all the cells that intersect the lobe
    std::vector<std::array<int, 2>> get_cells_intersecting_lobe( const Lobe & lobe );

    std::vector<std::pair<std::array<int, 2>, double>> compute_intersection( const Lobe & lobe, int N = 15 );

    // Figure out which pixel a given point is in, returning the indices of the lowest left corner
    std::array<int, 2> locate_point( const Vector2 & coordinates );
};

} // namespace Flowtastic