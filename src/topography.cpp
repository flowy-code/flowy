#include "topography.hpp"
#include "definitions.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xbuilder.hpp"
#include <algorithm>
#include <stdexcept>

namespace Flowtastic
{

std::array<int, 2> Topography::locate_point( const Vector2 & coordinates )
{
    const bool outside_x = coordinates[0] < x_data[0] || coordinates[0] >= x_data.periodic( -1 ) + cell_size();
    const bool outside_y = coordinates[1] < y_data[0] || coordinates[1] >= y_data.periodic( -1 ) + cell_size();

    if( outside_x || outside_y )
    {
        throw std::runtime_error( "Cannot locate point, because coordinates are outside of grid!" );
    }

    const int idx_x_lower = int( ( coordinates[0] - x_data[0] ) / cell_size() );
    const int idx_y_lower = int( ( coordinates[1] - y_data[0] ) / cell_size() );
    return { idx_x_lower, idx_y_lower };
}

bool Topography::point_in_cell( int idx_i, int idx_j, const Vector2 & point )
{
    // Get the center of the cell
    const Vector2 center_of_cell = { x_data[idx_i] + 0.5 * cell_size(), y_data[idx_j] + 0.5 * cell_size() };

    // Note that the left and bottom border of the cell are included, while the right and top are not
    // clang-format off
    return  (point[0] >= center_of_cell[0] -cell_size()/2.0) &&
            (point[0] < center_of_cell[0] + cell_size()/2.0) &&
            (point[1] >= center_of_cell[1] -cell_size()/2.0) &&
            (point[1] < center_of_cell[1] + cell_size()/2.0);
    // clang-format on
}

bool Topography::line_intersects_cell( int idx_i, int idx_j, double slope_xy, double offset )
{
    // Get the center of the cell
    const Vector2 center_of_cell = { x_data[idx_i] + 0.5 * cell_size(), y_data[idx_j] + 0.5 * cell_size() };

    // We have to check if the line intersects any ouf the four sides of the square

    auto check_x = [&]( double xc )
    { return xc >= center_of_cell[0] - 0.5 * cell_size() && xc < center_of_cell[0] + 0.5 * cell_size(); };

    auto check_y = [&]( double yc )
    { return yc >= center_of_cell[1] - 0.5 * cell_size() && yc < center_of_cell[1] + 0.5 * cell_size(); };

    double yc{}, xc{};

    // Check bottom
    yc = center_of_cell[1] - cell_size() / 2;
    xc = ( yc - offset ) / slope_xy;
    if( check_x( xc ) )
    {
        return true;
    }

    // Check top
    yc = center_of_cell[1] + cell_size() / 2;
    xc = ( yc - offset ) / slope_xy;
    if( check_x( xc ) )
    {
        return true;
    }

    // Check left
    xc = center_of_cell[0] - cell_size() / 2;
    yc = slope_xy * xc + offset;
    if( check_y( yc ) )
    {
        return true;
    }

    // Check right
    xc = center_of_cell[0] + cell_size() / 2;
    yc = slope_xy * xc + offset;
    return check_y( yc );
}

Topography::BoundingBox Topography::bounding_box( const Vector2 & center, double radius )
{
    const auto [idx_x_lower, idx_y_lower] = locate_point( center );
    int number_of_cells_to_include        = std::ceil( radius / cell_size() );

    Topography::BoundingBox res{};
    res.idx_x_lower  = std::clamp<int>( idx_x_lower - number_of_cells_to_include, 0, x_data.size() - 1 );
    res.idx_x_higher = std::clamp<int>( idx_x_lower + number_of_cells_to_include, 0, x_data.size() - 1 );
    res.idx_y_lower  = std::clamp<int>( idx_y_lower - number_of_cells_to_include, 0, y_data.size() - 1 );
    res.idx_y_higher = std::clamp<int>( idx_y_lower + number_of_cells_to_include, 0, y_data.size() - 1 );

    return res;
}

bool Topography::line_segment_intersects_cell( int idx_x, int idx_y, Vector2 x1, Vector2 x2 )
{

    const double slope  = ( x2[1] - x1[1] ) / ( x2[0] - x1[0] );
    const double offset = x1[1] - slope * x1[0];
    const double x_low  = std::min( x1[0], x2[0] );
    const double x_high = std::max( x1[0], x2[0] );

    if( x_data[idx_x] >= x_low && x_data[idx_x] < x_high )
    {
        return line_intersects_cell( idx_x, idx_y, slope, offset );
    }
    return false;
}

std::pair<MatrixX, Topography::BoundingBox> Topography::compute_intersection( const Lobe & lobe )
{
    auto bbox       = bounding_box( lobe.center, lobe.semi_axes[0] );
    constexpr int N = 15;

    const int n_rows = bbox.idx_y_higher - bbox.idx_y_lower + 1;
    const int n_cols = bbox.idx_x_higher - bbox.idx_x_lower + 1;

    MatrixX covered_fraction = xt::empty<double>( { n_rows, n_cols } );

    // We iterate over all cells, contained in the bounding box
    for( int idx_x = bbox.idx_x_lower; idx_x <= bbox.idx_x_higher; idx_x++ )
    {
        const double x_cell_low  = x_data[idx_x];
        const double x_cell_high = x_data[idx_x] + cell_size();
        const auto x_range_cell  = xt::linspace<double>( x_cell_low, x_cell_high, N );

        for( int idx_y = bbox.idx_y_lower; idx_y <= bbox.idx_y_higher; idx_y++ )
        {
            const double y_cell_low  = y_data[idx_y];
            const double y_cell_high = y_data[idx_y] + cell_size();
            const auto y_range_cell  = xt::linspace<double>( y_cell_low, y_cell_high, N );
            int n_hits               = 0;
            for( auto x : x_range_cell )
            {
                for( auto y : y_range_cell )
                {
                    if( lobe.is_point_in_lobe( { x, y } ) )
                    {
                        n_hits++;
                    }
                }
            }
            covered_fraction( idx_x - bbox.idx_x_lower, idx_y - bbox.idx_y_lower ) = double( n_hits ) / ( N * N );
        }
    }

    return { covered_fraction, bbox };
}

std::pair<double, Vector2> Topography::height_and_slope( const Vector2 & coordinates )
{
    const auto [idx_x_lower, idx_y_lower] = locate_point( coordinates );

    int idx_x_higher = idx_x_lower + 1;
    int idx_y_higher = idx_y_lower + 1;

    if( size_t( idx_x_higher ) == x_data.size() )
    {
        idx_x_higher = idx_x_lower;
    }

    if( size_t( idx_y_higher ) == y_data.size() )
    {
        idx_y_higher = idx_y_lower;
    }

    // We interpolate the height data with a function f(x,y) = alpha * x + beta * y + gamma * x * y + c
    // We want f(x,y) to coincide with the actual height data at the four enclosing corners
    // The height data is
    const double z00 = height_data( idx_x_lower, idx_y_lower );
    const double z01 = height_data( idx_x_higher, idx_y_lower );
    const double z10 = height_data( idx_x_lower, idx_y_higher );
    const double z11 = height_data( idx_x_higher, idx_y_higher );

    // And the values of x and y are
    const double x0 = x_data[idx_x_lower];
    const double x1 = x_data[idx_x_higher];
    const double y0 = x_data[idx_y_lower];
    const double y1 = x_data[idx_y_higher];

    // We transform the coordinates into the interval [0, 1]
    // with the transformations
    // x' -> (x - x0) / (x1 - x0)
    // y' -> (y - y0) / (y1 - y0)
    // In the new coordinates, this results in four equations,
    // from which the values of alpha, beta and gamma can be
    // derived:
    // f(0,0) = z00
    // f(1,0) = z10
    // f(0,1) = z01
    // f(1,1) = z11
    // c = f(0,0)
    // alpha = f(1,0) - f(0,0)
    // beta  = f(0,1) - f(0,0)
    // gamma = f(1,1) - f(1,0) + f(0,0) - f(0,1)
    const double c     = z00;
    const double alpha = z10 - z00;
    const double beta  = z01 - z00;
    const double gamma = z11 - z10 + z00 - z01;

    // Finally we get the desired value of z at the coordinatex (xc, yc)
    const double x_prime = ( coordinates[0] - x0 ) / ( x1 - x0 );
    const double y_prime = ( coordinates[1] - y0 ) / ( y1 - y0 );

    const double height = alpha * x_prime + beta * y_prime + gamma * x_prime * y_prime + c;

    // Now we compute the slope as
    // \nabla f(x,y) = [\del_x f(x,y), \del_y f(x,y)]
    // Derivative of x_prime wrt x
    const double x_prime_dx = 1.0 / ( x1 - x0 );
    // Derivative of y_prime wrt y
    const double y_prime_dy = 1.0 / ( y1 - y0 );

    const double slope_x = alpha * x_prime_dx + gamma * y_prime * x_prime_dx;
    const double slope_y = beta * y_prime_dy + gamma * x_prime * y_prime_dy;

    const Vector2 slope = { slope_x, slope_y };
    return { height, slope };
}

} // namespace Flowtastic