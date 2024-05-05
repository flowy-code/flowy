#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Flowy
{

struct LobeCells
{
    using cellvecT = std::vector<std::array<int, 2>>;
    cellvecT cells_intersecting{};
    cellvecT cells_enclosed{};

    static constexpr int n_trapz = 5;
    using trapzT                 = std::array<double, n_trapz>;
    std::vector<trapzT> cell_trapz_values{};
};

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

    Topography( const MatrixX & height_data, const VectorX & x_data, const VectorX & y_data, double no_data_value )
            : height_data( height_data ),
              hazard( xt::zeros_like( height_data ) ),
              x_data( x_data ),
              y_data( y_data ),
              no_data_value( no_data_value )
    {
        // Check that x_data spacing and y_data spacing is the same
        double delta_x = x_data[1] - x_data[0];
        double delta_y = y_data[1] - y_data[0];
        if( !xt::isclose( delta_x, delta_y )[0] )
        {
            throw std::runtime_error( fmt::format(
                "The spacing between x and y must be the same! It is delta_x = {} and delta_y = {}", delta_x,
                delta_y ) );
        }
    }

    Topography() = default;

    MatrixX height_data{}; // The heights of the cells
    MatrixX hazard{};      // Contains data on the cumulative descendents
    VectorX x_data{};
    VectorX y_data{};
    double no_data_value{};

    inline double get_height( int idx_x, int idx_y ) const
    {
        return height_data( idx_x, idx_y );
    }

    inline double get_height( const Vector2 & point ) const
    {
        auto [idx_x, idx_y] = locate_point( point );
        return height_data( idx_x, idx_y );
    }

    inline void set_height( int idx_x, int idx_y, double height )
    {
        height_data( idx_x, idx_y ) = height;
    }

    inline void set_height( const Vector2 & point, double height )
    {
        auto [idx_x, idx_y]         = locate_point( point );
        height_data( idx_x, idx_y ) = height;
    }

    inline double cell_size() const
    {
        return x_data[1] - x_data[0];
    };

    inline double volume() const
    {
        return xt::sum( height_data )() * cell_size() * cell_size();
    }

    inline double area( double thresh = 0.0 ) const
    {
        return xt::count_nonzero( xt::filter( height_data, height_data > thresh ) )() * cell_size() * cell_size();
    }

    // Calculate the height and the slope at coordinates
    // via linear interpolation from the square grid
    std::pair<double, Vector2> height_and_slope( const Vector2 & coordinates ) const;

    // Calculate the slope between points, given their heights
    double slope_between_points(
        const Vector2 & point1, const Vector2 & point2, std::optional<double> min_height_drop = 0.0 ) const;

    // Compute the indices of a rectangular bounding box
    // The box is computed such that a circle with centered at 'center' with radius 'radius'
    // Is completely contained in the bounding box
    BoundingBox bounding_box( const Vector2 & center, double extent_x, double extent_y ) const;

    // Find all the cells that intersect the lobe and all the cells that are fully enclosed by the lobe
    LobeCells get_cells_intersecting_lobe( const Lobe & lobe, std::optional<int> idx_cache = std::nullopt );

    double rasterize_cell_grid( int idx_x, int idx_y, const Lobe & lobe ) const;
    double rasterize_cell_trapz( LobeCells::trapzT & trapz_values ) const;

    // Find the fraction of the cells covered by the lobe by rasterizing each cell
    // into a grid of N*N points
    // This returns a vector of pairs
    // - the first entry of each pair contains an array<int, 2> with the idx_i, idx_j of the intersected cell
    // - the second entry contains the fraction of the cell that is covered by the ellips
    std::vector<std::pair<std::array<int, 2>, double>>
    compute_intersection( const Lobe & lobe, std::optional<int> idx_cache = std::nullopt );

    // Adds the lobe thickness to the topography, according to its fractional intersection with the cells
    void add_lobe( const Lobe & lobe, bool volume_correction, std::optional<int> idx_cache = std::nullopt );

    // Computes the hazard for a flow
    void compute_hazard_flow( const std::vector<Lobe> & lobes );

    // Check if a point is near the boundary
    bool is_point_near_boundary( const Vector2 & coordinates, double radius ) const;

    // Figure out which cell a given point is in, returning the indices of the lowest left corner
    std::array<int, 2> locate_point( const Vector2 & coordinates ) const;

    // Add a topography to the topography object
    // The filling_parameter is a scale factor to multiply the height by when adding to the topography
    // Apparently the filling_parameter accounts for "subsurface flows"
    void add_to_topography( const Topography & topography_to_a, double filling_parameter = 1.0 );

    Vector2 find_preliminary_budding_point( const Lobe & lobe, int npoints ) const;

    void reset_intersection_cache( int N );

private:
    std::vector<std::optional<LobeCells>> intersection_cache{};
};

} // namespace Flowy
