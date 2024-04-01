#include "definitions.hpp"
#include "lobe.hpp"
#include "math.hpp"
#include "topography.hpp"
#include "xtensor/xbuilder.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

TEST_CASE( "bounding_box", "[bounding_box]" )
{
    Flowtastic::VectorX x_data      = xt::arange<double>( 1.0, 20.0, 1.0 );
    Flowtastic::VectorX y_data      = xt::arange<double>( 5.0, 11.0, 1.0 );
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

TEST_CASE( "test_compute_intersection", "[intersection]" )
{
    Flowtastic::VectorX x_data      = xt::arange<double>( 0.0, 10.0, 1.0 );
    Flowtastic::VectorX y_data      = xt::arange<double>( 0.0, 10.0, 1.0 );
    Flowtastic::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowtastic::Topography( height_data, x_data, y_data );

    Flowtastic::Lobe my_lobe;
    my_lobe.center          = { 2.5, 5.5 };
    my_lobe.semi_axes       = { 2, 1 };
    my_lobe.azimuthal_angle = 0.0;

    auto [intersection, bbox] = topography.compute_intersection( my_lobe );

    fmt::print( "bbox.idx_x_higher = {}\n", bbox.idx_x_higher );
    fmt::print( "bbox.idx_x_lower = {}\n", bbox.idx_x_lower );
    fmt::print( "bbox.idx_y_lower = {}\n", bbox.idx_y_lower );
    fmt::print( "bbox.idx_y_higher = {}\n", bbox.idx_y_higher );
    fmt::print( "intersection = {}\n", fmt::streamed( intersection( 0, 0 ) ) );
}