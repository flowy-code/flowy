// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/topography.hpp"
#include "flowy/include/asc_file.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "thirdparty/tsl/robin_set.h"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>
#include <unordered_set>
#include <vector>

namespace Flowy
{

bool Topography::is_point_near_boundary( const Vector2 & coordinates, double radius ) const
{
    int n = std::ceil( radius / cell_size() );

    const bool near_x_boundary
        = coordinates[0] < x_data[0] + n * cell_size() || coordinates[0] >= x_data.periodic( -1 ) - n * cell_size();
    const bool near_y_boundary
        = coordinates[1] < y_data[0] + n * cell_size() || coordinates[1] >= y_data.periodic( -1 ) - n * cell_size();
    return near_x_boundary || near_y_boundary;
}

std::array<int, 2> Topography::locate_point( const Vector2 & coordinates ) const
{
    const bool outside_x = coordinates[0] < x_data[0] || coordinates[0] >= x_data.periodic( -1 ) + cell_size();
    const bool outside_y = coordinates[1] < y_data[0] || coordinates[1] >= y_data.periodic( -1 ) + cell_size();

    if( outside_x || outside_y )
    {
        throw std::runtime_error( "Cannot locate point, because coordinates are outside of grid!" );
    }

    const int idx_x = static_cast<int>( ( coordinates[0] - x_data[0] ) / cell_size() );
    const int idx_y = static_cast<int>( ( coordinates[1] - y_data[0] ) / cell_size() );
    return { idx_x, idx_y };
}

Topography::BoundingBox Topography::bounding_box( const Vector2 & center, double extent_x, double extent_y ) const
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

// Helper struct
struct RowIntersectionData
{
    int idx_row{};

    double y_top;
    double y_bot;

    std::optional<double> x_left_top{};
    std::optional<double> x_right_top{};
    std::optional<double> x_left_bot{};
    std::optional<double> x_right_bot{};

    std::optional<int> idx_x_left_top{};
    std::optional<int> idx_x_right_top{};
    std::optional<int> idx_x_left_bot{};
    std::optional<int> idx_x_right_bot{};

    // vectors of intermediate x and y values (used for trapezoidal rule)
    std::vector<double> x_left_arr{};
    std::vector<double> x_right_arr{};
    std::vector<double> y_left_arr{};
};

LobeCells Topography::get_cells_intersecting_lobe( const Lobe & lobe, std::optional<int> idx_cache )
{
    // Can we use the cache?
    bool use_cache
        = idx_cache.has_value() && ( intersection_cache.size() > static_cast<std::size_t>( idx_cache.value() ) );

    if( use_cache )
    {
        std::optional<LobeCells> cached_lobe_cells = intersection_cache[idx_cache.value()];

        // Does the cache already contain a value?
        if( cached_lobe_cells.has_value() ) // If yes, we return it
        {
            return cached_lobe_cells.value();
        }
    }

    LobeCells res{};

    auto extent_xy = lobe.extent_xy();

    int idx_y_min    = ( lobe.center[1] - extent_xy[1] - y_data[0] ) / cell_size();
    int idx_y_max    = ( lobe.center[1] + extent_xy[1] - y_data[0] ) / cell_size();
    const int n_rows = idx_y_max - idx_y_min + 1;

    auto row_data = std::vector<RowIntersectionData>( n_rows );

    // Measure the row data. We can skip the first row, since we know there are no intersections there
    for( int irow = 0; irow < n_rows; irow++ )
    {
        // y-value at the bottom of the row
        const double y = y_data[idx_y_min + irow];

        row_data[irow].idx_row = irow;
        row_data[irow].y_bot   = y;
        row_data[irow].y_top   = y_data[idx_y_min + irow + 1];

        if( irow == 0 )
            continue;

        // These two points define a horizontal line, at the bottom of the current row
        const Vector2 x1  = { lobe.center[0] - extent_xy[0], y };
        const Vector2 x2  = { lobe.center[0] + extent_xy[0], y };
        const auto points = lobe.line_segment_intersects( x1, x2 );

        // Only if intersections are found, we can unpack them into x indices
        if( points.has_value() )
        {
            const auto p1 = points.value()[0];
            const auto p2 = points.value()[1];

            // Rows share a bottom and a top, so we only need one linesegment intersects
            row_data[irow - 1].x_left_top  = p1[0];
            row_data[irow - 1].x_right_top = p2[0];
            row_data[irow].x_left_bot      = p1[0];
            row_data[irow].x_right_bot     = p2[0];

            row_data[irow - 1].idx_x_left_top  = ( p1[0] - x_data[0] ) / cell_size();
            row_data[irow - 1].idx_x_right_top = ( p2[0] - x_data[0] ) / cell_size();

            row_data[irow].idx_x_left_bot  = ( p1[0] - x_data[0] ) / cell_size();
            row_data[irow].idx_x_right_bot = ( p2[0] - x_data[0] ) / cell_size();
        }
    }

    // push_back enclosed cells with x index in the interval [idx_start, idx_stop]
    auto push_back_enclosed_cells = [&]( int idx_start, int idx_stop, int idx_y )
    {
        for( int idx_x = idx_start; idx_x <= idx_stop; idx_x++ )
        {
            res.cells_enclosed.push_back( { idx_x, idx_y } );
        }
    };

    // push_back intersected cells with x index in the interval [idx_start, idx_stop]
    auto push_back_intersected_cells
        = [&]( int idx_start, int idx_stop, int idx_y,
               const std::array<std::optional<std::array<double, 2>>, LobeCells::n_trapz> & intersection_values )
    {
        for( int idx_x = idx_start; idx_x <= idx_stop; idx_x++ )
        {
            res.cells_intersecting.push_back( { idx_x, idx_y } );
            LobeCells::trapzT trapz_values{};

            double x_cell_min = x_data[idx_x];
            double x_cell_max = x_data[idx_x] + cell_size();

            // Figure out the trapz value
            for( size_t i = 0; i < intersection_values.size(); i++ )
            {
                const auto & inter = intersection_values[i];

                if( inter.has_value() )
                {
                    const double x_left  = std::clamp( inter.value()[0], x_cell_min, x_cell_max );
                    const double x_right = std::clamp( inter.value()[1], x_cell_min, x_cell_max );
                    trapz_values[i]      = x_right - x_left;
                }
                else
                {
                    trapz_values[i] = 0;
                }
            }
            res.cell_trapz_values.push_back( trapz_values );
        }
    };

    for( int irow = 0; irow < n_rows; irow++ )
    {
        const int idx_y = idx_y_min + irow;

        const auto & row_data_cur = row_data[irow];

        std::array<std::optional<std::array<double, 2>>, LobeCells::n_trapz> intersections{};

        // Fill in the already known intersections (at the bottom and the top of the row)
        if( row_data_cur.x_left_bot.has_value() && row_data_cur.x_right_bot.has_value() )
            intersections[0] = { row_data_cur.x_left_bot.value(), row_data_cur.x_right_bot.value() };

        if( row_data_cur.x_left_top.has_value() && row_data_cur.x_right_top.has_value() )
            intersections[LobeCells::n_trapz - 1]
                = { row_data_cur.x_left_top.value(), row_data_cur.x_right_top.value() };

        // Loop goes from i=1 to i=N-2, since the top and the bottom of the row are already known
        for( int i = 1; i < LobeCells::n_trapz - 1; i++ )
        {
            const double y_cur = y_data[idx_y] + cell_size() / ( LobeCells::n_trapz - 1 ) * i;

            // Get the line segment intersections
            const Vector2 x1  = { lobe.center[0] - extent_xy[0], y_cur };
            const Vector2 x2  = { lobe.center[0] + extent_xy[0], y_cur };
            const auto points = lobe.line_segment_intersects( x1, x2 );

            if( points.has_value() )
            {
                intersections[i] = { points.value()[0][0], points.value()[1][0] };
            }
        }

        // We treat the first and the last row separately, since here, there are no intersections
        // with the previous row (in case of the first) or the current row (in case of the last row)
        if( irow == 0 )
        {
            push_back_intersected_cells(
                row_data_cur.idx_x_left_top.value(), row_data_cur.idx_x_right_top.value(), idx_y_min, intersections );
        }
        else if( irow == n_rows - 1 )
        {
            push_back_intersected_cells(
                row_data_cur.idx_x_left_bot.value(), row_data_cur.idx_x_right_bot.value(), idx_y_max, intersections );
        }
        else
        {
            const int start_left
                = std::min<int>( row_data_cur.idx_x_left_bot.value(), row_data_cur.idx_x_left_top.value() );
            const int stop_left
                = std::max<int>( row_data_cur.idx_x_left_bot.value(), row_data_cur.idx_x_left_top.value() );
            push_back_intersected_cells( start_left, stop_left, idx_y, intersections );

            int start_right
                = std::min<int>( row_data_cur.idx_x_right_bot.value(), row_data_cur.idx_x_right_top.value() );
            const int stop_right
                = std::max<int>( row_data_cur.idx_x_right_bot.value(), row_data_cur.idx_x_right_top.value() );

            // If stop_left and start_right co-incide, which can happen when the tip of an ellipse barely touches a row
            // We need to make sure not to double count an intersected cell
            if( stop_left == start_right )
                start_right++;

            push_back_intersected_cells( start_right, stop_right, idx_y, intersections );
            push_back_enclosed_cells( stop_left + 1, start_right - 1, idx_y );
        }
    }

    // If the cache is used, we copy the intersection data there
    if( use_cache )
    {
        intersection_cache[idx_cache.value()] = res;
    }

    return res;
}

double Topography::rasterize_cell_grid( int idx_x, int idx_y, const Lobe & lobe ) const
{
    constexpr int N = 15;

    const double cell_size = this->cell_size();
    const double step      = cell_size / ( N - 1 );

    int n_hits = 0;
    for( int ix = 0; ix < N; ix++ )
    {
        const double x = x_data[idx_x] + step * ix;
        for( int iy = 0; iy < N; iy++ )
        {
            const double y = y_data[idx_y] + step * iy;
            if( lobe.is_point_in_lobe( { x, y } ) )
            {
                n_hits++;
            }
        }
    }

    const double fraction = static_cast<double>( n_hits ) / static_cast<double>( N * N );
    return fraction;
}

double Topography::rasterize_cell_trapz( LobeCells::trapzT & trapz_values ) const
{
    const double cell_size     = this->cell_size();
    const double cell_area     = cell_size * cell_size;
    const double trapz_spacing = cell_size / ( LobeCells::n_trapz - 1 );

    // Evaluate trapezoidal rule
    double area = 0;
    for( size_t ix = 0; ix < trapz_values.size() - 1; ix++ )
    {
        area += 0.5 * ( trapz_values[ix] + trapz_values[ix + 1] );
    }

    const double fraction = area * trapz_spacing / cell_area;
    return fraction;
}

std::vector<std::pair<std::array<int, 2>, double>>
Topography::compute_intersection( const Lobe & lobe, std::optional<int> idx_cache )
{
    auto lobe_cells = get_cells_intersecting_lobe( lobe, idx_cache );

    std::vector<std::pair<std::array<int, 2>, double>> res{};
    res.reserve( lobe_cells.cells_intersecting.size() + lobe_cells.cells_enclosed.size() );

    // All enclosed cells are fully covered
    for( const auto & [idx_x, idx_y] : lobe_cells.cells_enclosed )
    {
        res.push_back( { { idx_x, idx_y }, 1.0 } );
    }

    // The intersecting cells get rasterized
    for( size_t i = 0; i < lobe_cells.cells_intersecting.size(); i++ )
    {
        const auto & [idx_x, idx_y] = lobe_cells.cells_intersecting[i];
        const double fraction       = rasterize_cell_trapz( lobe_cells.cell_trapz_values[i] );
        res.push_back( { { idx_x, idx_y }, fraction } );
    }

    return res;
}

struct hash_pair
{
    template<class T>
    size_t operator()( const std::array<T, 2> & p ) const
    {
        size_t seed = 0;
        for( T i : p )
        {
            seed ^= std::hash<T>{}( i ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
        }
        return seed;
    }
};

void Topography::compute_hazard_flow( const std::vector<Lobe> & lobes )
{
    tsl::robin_set<std::array<int, 2>, hash_pair> parent_set{};

    // This computes the hazard for *one* flow
    // For one flow, the hazard of a cell is the maximum of lobe.n_descendant over all lobes touching it
    for( size_t idx = 0; idx < lobes.size(); idx++ )
    {
        const auto & lobe     = lobes[idx];
        const auto lobe_cells = get_cells_intersecting_lobe( lobe, idx );

        parent_set.clear();

        if( lobe.idx_parent.has_value() )
        {
            const auto & lobe_parent     = lobes[lobe.idx_parent.value()];
            const auto lobe_cells_parent = get_cells_intersecting_lobe( lobe_parent, lobe.idx_parent.value() );

            for( const auto & cells : { lobe_cells_parent.cells_enclosed, lobe_cells_parent.cells_intersecting } )
            {
                for( const auto & c : cells )
                {
                    parent_set.insert( c );
                }
            }
        }

        for( const auto & cells : { lobe_cells.cells_enclosed, lobe_cells.cells_intersecting } )
        {
            for( const auto & [idx_x, idx_y] : cells )
            {
                if( parent_set.contains( { idx_x, idx_y } ) )
                    continue;
                hazard( idx_x, idx_y ) += lobe.n_descendents + 1;
            }
        }
    }
}

std::pair<double, Vector2> Topography::height_and_slope( const Vector2 & coordinates ) const
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

    if( Z00 == no_data_value || Z10 == no_data_value || Z01 == no_data_value || Z11 == no_data_value )
    {
        return { no_data_value, { no_data_value, no_data_value } };
    }

    const double alpha = Z10 - Z00;
    const double beta  = Z01 - Z00;
    const double gamma = Z11 + Z00 - Z10 - Z01;

    const Vector2 xp = ( coordinates - cell_center_lower_left ) / cell_size();

    const double height = Z00 + alpha * xp[0] + beta * xp[1] + gamma * xp[0] * xp[1];
    const Vector2 slope = { alpha + gamma * xp[1], beta + gamma * xp[0] };

    return { height, -slope / cell_size() };
}

double Topography::slope_between_points(
    const Vector2 & point1, const Vector2 & point2, std::optional<double> min_height_drop ) const
{
    const double height1 = height_and_slope( point1 ).first;
    const double height2 = height_and_slope( point2 ).first;

    const Vector2 diff       = point2 - point1;
    double norm              = std::sqrt( diff[0] * diff[0] + diff[1] * diff[1] );
    double height_difference = -( height2 - height1 );

    if( min_height_drop.has_value() )
    {
        height_difference = std::max( min_height_drop.value(), height_difference );
    }

    double slope = height_difference / norm;
    return slope;
}

void Topography::add_lobe( const Lobe & lobe, bool volume_correction, std::optional<int> idx_cache )
{
    // In this function we simply add the thickness of the lobe to the topography
    // First, we find the intersected cells and the covered fractions

    std::vector<std::pair<std::array<int, 2>, double>> intersection_data = compute_intersection( lobe, idx_cache );
    double volume_added            = 0.0; // Volume added to the topography from rasterization
    double area_intersecting_cells = 0.0; // Total area covered by intersecting cells
    const double cell_area         = cell_size() * cell_size();

    // Then we add the tickness according to the fractions
    for( auto const & [indices, fraction] : intersection_data )
    {
        const double cell_height = fraction * lobe.thickness;
        height_data( indices[0], indices[1] ) += cell_height;
        volume_added += cell_height * cell_area;
        if( fraction < 1.0 )
        {
            area_intersecting_cells += fraction * cell_area;
        }
    }

    // Optionally compute the volume correction
    if( volume_correction )
    {
        const double volume_to_add     = lobe.volume() - volume_added;
        const double avg_height_to_add = volume_to_add / area_intersecting_cells;
        // Go over the cells
        for( auto const & [indices, fraction] : intersection_data )
        {
            if( fraction < 1.0 )
            {
                const double cell_height = fraction * avg_height_to_add;
                height_data( indices[0], indices[1] ) += cell_height;
            }
        }
    }
}

Vector2 Topography::find_preliminary_budding_point( const Lobe & lobe, int npoints ) const
{
    // First, we rasterize the perimeter of the ellipse
    std::vector<Vector2> perimeter = lobe.rasterize_perimeter( npoints );

    // Then, we find the point of minimal elevation amongst the rasterized points on the perimeter
    auto min_elevation_point_it = std::min_element(
        perimeter.begin(), perimeter.end(), [&]( const Vector2 & p1, const Vector2 & p2 )
        { return height_and_slope( p1 ).first < height_and_slope( p2 ).first; } );

    return *min_elevation_point_it;
}

void Topography::reset_intersection_cache( int N )
{
    intersection_cache = std::vector<std::optional<LobeCells>>( N, std::nullopt );
}

void Topography::add_to_topography( const Topography & topography_to_add, double filling_parameter )
{
    // Loop over the cells in the topography
    for( size_t idx_x = 0; idx_x < this->x_data.size(); idx_x++ )
    {
        for( size_t idx_y = 0; idx_y < this->y_data.size(); idx_y++ )
        {
            // Get the point
            const Vector2 point = { this->x_data[idx_x], this->y_data[idx_y] };

            // Skip if this point is outside the extents of topography_to_add
            if( topography_to_add.is_point_near_boundary( point, 0.0 ) )
            {
                continue;
            }

            auto [height_to_add, slope] = topography_to_add.height_and_slope( point );
            auto old_height             = height_data( idx_x, idx_y );

            if( height_to_add == topography_to_add.no_data_value || old_height == this->no_data_value )
                continue;

            // Add this height to the current topography
            height_data( idx_x, idx_y ) += filling_parameter * height_to_add;
        }
    }
}

} // namespace Flowy
