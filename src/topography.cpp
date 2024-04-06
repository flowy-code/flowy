#include "topography.hpp"
#include "asc_file.hpp"
#include "definitions.hpp"
#include "xtensor/xbuilder.hpp"
#include <fmt/ranges.h>
#include <algorithm>
#include <stdexcept>

namespace Flowy
{

AscFile Topography::to_asc_file()
{
    AscFile asc_file{};
    asc_file.lower_left_corner = { x_data[0], y_data[0] };
    asc_file.cell_size         = cell_size();
    asc_file.height_data       = height_data;
    return asc_file;
}

bool Topography::is_point_near_boundary( const Vector2 & coordinates, double radius )
{
    int n = std::ceil( radius / cell_size() );

    const bool near_x_boundary
        = coordinates[0] < x_data[0] + n * cell_size() || coordinates[0] >= x_data.periodic( -1 ) - n * cell_size();
    const bool near_y_boundary
        = coordinates[1] < y_data[0] + n * cell_size() || coordinates[1] >= y_data.periodic( -1 ) - n * cell_size();
    return near_x_boundary || near_y_boundary;
}

std::array<int, 2> Topography::locate_point( const Vector2 & coordinates )
{
    const bool outside_x = coordinates[0] < x_data[0] || coordinates[0] >= x_data.periodic( -1 ) + cell_size();
    const bool outside_y = coordinates[1] < y_data[0] || coordinates[1] >= y_data.periodic( -1 ) + cell_size();

    if( outside_x || outside_y )
    {
        throw std::runtime_error( "Cannot locate point, because coordinates are outside of grid!" );
    }

    const int idx_x = int( ( coordinates[0] - x_data[0] ) / cell_size() );
    const int idx_y = int( ( coordinates[1] - y_data[0] ) / cell_size() );
    return { idx_x, idx_y };
}

Topography::BoundingBox Topography::bounding_box( const Vector2 & center, double extent_x, double extent_y )
{
    const auto [idx_x_lower, idx_y_lower]  = locate_point( center );
    const int number_of_cells_to_include_x = std::ceil( extent_x / cell_size() );
    const int number_of_cells_to_include_y = std::ceil( extent_y / cell_size() );

    Topography::BoundingBox res{};
    res.idx_x_lower  = std::clamp<int>( idx_x_lower - number_of_cells_to_include_x, 0, x_data.size() - 1 );
    res.idx_x_higher = std::clamp<int>( idx_x_lower + number_of_cells_to_include_x, 0, x_data.size() - 1 );
    res.idx_y_lower  = std::clamp<int>( idx_y_lower - number_of_cells_to_include_y, 0, y_data.size() - 1 );
    res.idx_y_higher = std::clamp<int>( idx_y_lower + number_of_cells_to_include_y, 0, y_data.size() - 1 );

    return res;
}

std::pair<std::vector<std::array<int, 2>>, std::vector<std::array<int, 2>>>
Topography::get_cells_intersecting_lobe( const Lobe & lobe )
{
    std::vector<std::array<int, 2>> cells_intersecting{};
    std::vector<std::array<int, 2>> cells_enclosed{};

    // First, we find all candidates with the bounding_box function from the topography
    auto extent_xy  = lobe.extent_xy();
    const auto bbox = bounding_box( lobe.center, extent_xy[0], extent_xy[1] );

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
                   lobe.line_segment_intersects( point_lb,point_rb )
                || lobe.line_segment_intersects( point_rb,point_rt )
                || lobe.line_segment_intersects( point_rt,point_lt )
                || lobe.line_segment_intersects( point_lt,point_lb )
            )
            {
                cells_intersecting.push_back( { idx_x, idx_y } );
            } else {
                if (lobe.is_point_in_lobe(point_lb))
                    cells_enclosed.push_back( {idx_x, idx_y} );
            }
            // clang-format on
        }
    }

    return { cells_intersecting, cells_enclosed };
}

std::vector<std::pair<std::array<int, 2>, double>> Topography::compute_intersection( const Lobe & lobe, int N )
{
    const double N2                           = N * N;
    auto [cells_intersecting, cells_enclosed] = get_cells_intersecting_lobe( lobe );

    std::vector<std::pair<std::array<int, 2>, double>> res{};
    res.reserve( cells_intersecting.size() + cells_enclosed.size() );

    const double step = cell_size() / N;

    for( const auto & [idx_x, idx_y] : cells_enclosed )
    {
        res.push_back( { { idx_x, idx_y }, 1.0 } );
    }

    for( const auto & [idx_x, idx_y] : cells_intersecting )
    {
        double n_hits = 0.0;

        for( int ix = 0; ix < N; ix++ )
        {
            const double x = x_data[idx_x] + step * ix;
            for( int iy = 0; iy < N; iy++ )
            {
                const double y = y_data[idx_y] + step * iy;
                if( lobe.is_point_in_lobe( { x, y } ) )
                {
                    n_hits += 1.0;
                }
            }
        }
        const double fraction = n_hits / N2;
        res.push_back( { { idx_x, idx_y }, fraction } );
    }

    return res;
}

std::pair<double, Vector2> Topography::height_and_slope( const Vector2 & coordinates )
{
    const auto [idx_x, idx_y] = locate_point( coordinates );
    const Vector2 cell_center = { x_data[idx_x] + 0.5 * cell_size(), y_data[idx_y] + 0.5 * cell_size() };

    int idx_x_lower{}, idx_x_higher{};
    int idx_y_lower{}, idx_y_higher{};

    if( coordinates[0] > cell_center[0] )
    {
        idx_x_lower  = idx_x;
        idx_x_higher = std::min<int>( idx_x + 1, x_data.size() - 1 );
    }
    else
    {
        idx_x_lower  = std::max<int>( idx_x - 1, 0 );
        idx_x_higher = idx_x;
    }

    if( coordinates[1] > cell_center[1] )
    {
        idx_y_lower  = idx_y;
        idx_y_higher = std::min<int>( idx_y + 1, y_data.size() - 1 );
    }
    else
    {
        idx_y_lower  = std::max<int>( idx_y - 1, 0 );
        idx_y_higher = idx_y;
    }

    const Vector2 cell_center_lower_left
        = { x_data[idx_x_lower] + 0.5 * cell_size(), y_data[idx_y_lower] + 0.5 * cell_size() };

    const double Z00 = height_data( idx_x_lower, idx_y_lower );
    const double Z10 = height_data( idx_x_higher, idx_y_lower );
    const double Z01 = height_data( idx_x_lower, idx_y_higher );
    const double Z11 = height_data( idx_x_higher, idx_y_higher );

    const double alpha = Z10 - Z00;
    const double beta  = Z01 - Z00;
    const double gamma = Z11 + Z00 - Z10 - Z01;

    const Vector2 xp = ( coordinates - cell_center_lower_left ) / cell_size();

    const double height = Z00 + alpha * xp[0] + beta * xp[1] + gamma * xp[0] * xp[1];
    const Vector2 slope = { alpha + gamma * xp[1], beta + gamma * xp[0] };

    return { height, slope/cell_size() };
}

void Topography::add_lobe( const Lobe & lobe )
{
    // In this function we simply add the thickness of the lobe to the topography
    // First, we find the intersected cells and the covered fractions
    std::vector<std::pair<std::array<int, 2>, double>> intersection_data = compute_intersection( lobe );

    // Then we add the tickness according to the fractions
    for( auto const & [indices, fraction] : intersection_data )
    {
        height_data( indices[0], indices[1] ) += fraction * lobe.thickness;
    }
}

Vector2 Topography::find_preliminary_budding_point( const Lobe & lobe, int npoints )
{
    // First, we rasterize the perimeter of the ellipse
    std::vector<Vector2> perimeter = lobe.rasterize_perimeter( npoints );

    // Then, we find the point of minimal elevation amongst the rasterized points on the perimeter
    auto min_elevation_point_it = std::min_element(
        perimeter.begin(), perimeter.end(), [&]( const Vector2 & p1, const Vector2 & p2 )
        { return height_and_slope( p1 ).first < height_and_slope( p2 ).first; } );

    return *min_elevation_point_it;
}

} // namespace Flowy