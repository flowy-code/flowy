#include "topography.hpp"
#include "definitions.hpp"
#include "xtensor/xbuilder.hpp"
#include <fmt/ranges.h>
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

Topography::BoundingBox Topography::bounding_box( const Vector2 & center, double radius )
{
    const auto [idx_x_lower, idx_y_lower] = locate_point( center );
    const int number_of_cells_to_include  = std::ceil( radius / cell_size() );

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
                   lobe.line_segment_intersects( point_lb,point_rb )
                || lobe.line_segment_intersects( point_rb,point_rt )
                || lobe.line_segment_intersects( point_rt,point_lt )
                || lobe.line_segment_intersects( point_lt,point_lb )
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
    const double N2         = N * N;
    auto cells_to_rasterize = get_cells_intersecting_lobe( lobe );

    std::vector<std::pair<std::array<int, 2>, double>> res{};
    res.reserve( cells_to_rasterize.size() );

    const double step = cell_size() / N;

    for( const auto & [idx_x, idx_y] : cells_to_rasterize )
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

    // We use central finite difference to compute the slope at the cell
    int idx_x_lower  = idx_x - 1;
    int idx_y_lower  = idx_y - 1;
    int idx_x_higher = idx_x + 1;
    int idx_y_higher = idx_y + 1;

    if( idx_x_lower < 0 )
    {
        idx_x_lower = 0;
    }

    if( idx_y_lower < 0 )
    {
        idx_y_lower = 0;
    }

    if( size_t( idx_x_higher ) == x_data.size() )
    {
        idx_x_higher = idx_x;
    }

    if( size_t( idx_y_higher ) == y_data.size() )
    {
        idx_y_higher = idx_y;
    }

    Vector2 slope{};
    if( coordinates[0] > cell_center[0] )
        slope[0] = ( height_data( idx_x_higher, idx_y ) - height_data( idx_x, idx_y ) ) / cell_size();
    else
        slope[0] = ( height_data( idx_x, idx_y ) - height_data( idx_x_lower, idx_y ) ) / cell_size();

    if( coordinates[1] > cell_center[1] )
        slope[1] = ( height_data( idx_x, idx_y_higher ) - height_data( idx_x, idx_y ) ) / cell_size();
    else
        slope[1] = ( height_data( idx_x, idx_y ) - height_data( idx_x, idx_y_lower ) ) / cell_size();

    // contribution of slope to height
    const Vector2 diff              = coordinates - cell_center;
    const double slope_contribution = slope[0] * diff[0] + slope[1] * diff[1];

    // For the height we linearly interpolate
    const double height = height_data( idx_x, idx_y ) + slope_contribution;

    return { height, slope };
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

} // namespace Flowtastic