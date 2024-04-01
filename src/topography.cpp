#include "topography.hpp"
#include "definitions.hpp"
#include "fmt/ostream.h"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xiterator.hpp"
#include <fmt/ranges.h>
#include <algorithm>
#include <limits>
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

std::vector<std::array<int, 2>> Topography::get_cells_intersecting_lobe( const Lobe & lobe )
{
    std::vector<std::array<int, 2>> res{};

    // First, we find all candidates with the bounding_box function from the topography
    const auto bbox = bounding_box( lobe.center, lobe.semi_axes[0] );

    // Now we test all the candidate cells from the big bounding box for intersection with the lobe
    for( int idx_x = bbox.idx_x_lower; idx_x <= bbox.idx_x_higher; idx_x++ )
    {
        for( int idx_y = bbox.idx_y_lower; idx_y <= bbox.idx_y_higher; idx_y++ )
        {
            // We check all corners of each cell
            Vector2 point_lb = { x_data[idx_x], y_data[idx_y] };
            Vector2 point_lt = { x_data[idx_x], y_data[idx_y] + cell_size() };
            Vector2 point_rb = { x_data[idx_x] + cell_size(), y_data[idx_y] };
            Vector2 point_rt = { x_data[idx_x] + cell_size(), y_data[idx_y] + cell_size() };
            // clang-format off
            if(
                   lobe.is_point_in_lobe( point_lb )
                || lobe.is_point_in_lobe( point_lt )
                || lobe.is_point_in_lobe( point_rb )
                || lobe.is_point_in_lobe( point_rt )
            )
            {
                res.push_back( { idx_x, idx_y } );
            }
            // clang-format on
        }
    }

    return res;
}

std::vector<std::pair<std::array<int, 2>, double>> Topography::compute_intersection( const Lobe & lobe, int N )
{
    std::vector<std::pair<std::array<int, 2>, double>> res{};

    auto cells_to_rasterize = get_cells_intersecting_lobe( lobe );

    for( auto [idx_x, idx_y] : cells_to_rasterize )
    {
        const auto x_range_cell = xt::linspace<double>( x_data[idx_x], x_data[idx_x] + cell_size(), N );
        const auto y_range_cell = xt::linspace<double>( y_data[idx_y], y_data[idx_y] + cell_size(), N );
        int n_hits              = 0;
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
        double fraction = double( n_hits ) / double( N * N );
        res.push_back( { { idx_x, idx_y }, fraction } );
    }
    return res;
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