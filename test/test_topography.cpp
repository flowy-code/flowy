#include "definitions.hpp"
#include "fmt/core.h"
#include "lobe.hpp"
#include "math.hpp"
#include "topography.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xio.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <vector>

TEST_CASE( "bounding_box", "[bounding_box]" )
{
    Flowtastic::VectorX x_data      = xt::arange<double>( 0.5, 19.5, 1.0 );
    Flowtastic::VectorX y_data      = xt::arange<double>( 4.5, 10.5, 1.0 );
    Flowtastic::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowtastic::Topography( height_data, x_data, y_data );

    Flowtastic::Vector2 point = { 11.4, 7.6 };

    int idx_point_x_expected  = 10;
    int idx_point_y_expected  = 3;
    double cell_size_expected = 1.0;

    auto point_indices = topography.locate_point( point );
    REQUIRE( idx_point_x_expected == point_indices[0] );
    REQUIRE( idx_point_y_expected == point_indices[1] );
    REQUIRE_THAT( topography.cell_size(), Catch::Matchers::WithinRel( cell_size_expected ) );

    fmt::print( "point.idx_x = {}, point.idx_y = {}\n", point_indices[0], point_indices[1] );

    // Note: this bounding box is chosen such that no clamping will occur
    auto bbox = topography.bounding_box( point, 2 );
    fmt::print( "bbox.idx_x_higher = {}\n", bbox.idx_x_higher );
    fmt::print( "bbox.idx_x_lower = {}\n", bbox.idx_x_lower );
    fmt::print( "bbox.idx_y_lower = {}\n", bbox.idx_y_lower );
    fmt::print( "bbox.idx_y_higher = {}\n", bbox.idx_y_higher );

    int idx_x_lower_expected  = idx_point_x_expected - 2;
    int idx_x_higher_expected = idx_point_x_expected + 2;
    int idx_y_lower_expected  = idx_point_y_expected - 2;
    int idx_y_higher_expected = idx_point_y_expected + 2;

    REQUIRE( bbox.idx_x_higher == idx_x_higher_expected );
    REQUIRE( bbox.idx_x_lower == idx_x_lower_expected );
    REQUIRE( bbox.idx_y_lower == idx_y_lower_expected );
    REQUIRE( bbox.idx_y_higher == idx_y_higher_expected );
}

TEST_CASE( "get_cells_intersecting_lobe", "[intersecting_lobe_cells]" )
{
    Flowtastic::VectorX x_data      = xt::arange<double>( -2.0, 2.0, 1.0 );
    Flowtastic::VectorX y_data      = xt::arange<double>( -2.0, 2.0, 1.0 );
    Flowtastic::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } ); // not really needed here

    auto topography = Flowtastic::Topography( height_data, x_data, y_data );

    Flowtastic::Lobe my_lobe;
    my_lobe.center          = { 0.0, 0.0 };
    my_lobe.semi_axes       = { 1.95, 0.95 };
    my_lobe.azimuthal_angle = Flowtastic::Math::pi / 2.0;

    // clang-format off
    std::vector<std::array<int, 2>> cell_indices_expected = { 
        { 1, 0 }, 
        { 1, 1 }, 
        { 1, 2 },
        { 1, 3 },
        { 2, 0 }, 
        { 2, 1 }, 
        { 2, 2 }, 
        { 2, 3 } };
    // clang-format on

    auto cell_indices = topography.get_cells_intersecting_lobe( my_lobe );

    for( size_t i = 0; i < cell_indices.size(); i++ )
    {
        fmt::print( " {} \n", cell_indices[i] );
        REQUIRE_THAT( cell_indices[i], Catch::Matchers::RangeEquals( cell_indices_expected[i] ) );
    }
}

TEST_CASE( "test_compute_intersection", "[intersection]" )
{
    Flowtastic::VectorX x_data      = xt::arange<double>( -3, 3, 1.0 );
    Flowtastic::VectorX y_data      = xt::arange<double>( -3, 3, 1.0 );
    Flowtastic::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowtastic::Topography( height_data, x_data, y_data );

    Flowtastic::Lobe my_lobe;
    my_lobe.center          = { 0, 0 };
    my_lobe.semi_axes       = { 1 - 1e-14, 1 - 1e-14 };
    my_lobe.azimuthal_angle = 0.0;

    // clang-format off
    std::vector<std::array<int, 2>> cell_indices_expected = { 
        { 2, 2 }, 
        { 2, 3 }, 
        { 3, 2 },
        { 3, 3 } };
    // clang-format on
    double expected_area_fraction = Flowtastic::Math::pi / 4.0;

    auto intersection_data = topography.compute_intersection( my_lobe, 30 );

    for( size_t i = 0; i < cell_indices_expected.size(); i++ )
    {
        auto indices  = intersection_data[i].first;
        auto fraction = intersection_data[i].second;

        fmt::print( " indices = {}, fraction = {}, fraction_expected {}\n", indices, fraction, expected_area_fraction );
        REQUIRE_THAT( indices, Catch::Matchers::RangeEquals( cell_indices_expected[i] ) );
        REQUIRE_THAT( fraction, Catch::Matchers::WithinRel( expected_area_fraction, 5e-2 ) );
    }
}