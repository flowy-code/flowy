// GPL v3 License
// Copyright 2023--present Flowy developers
#include "catch2/matchers/catch_matchers.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/topography.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <xtensor/xio.hpp>

TEST_CASE( "height_and_slope_test", "[topography]" )
{
    auto topography = Flowy::Topography();
    // x and y axes have the usual meaning: we assume that the y axis is already "flipped" from the ASC file default
    // (top to down)
    topography.height_data = { { 4.0, 1.0 }, { 4.0, 1.0 } };

    topography.x_data = { 1.0, 2.0 };
    topography.y_data = { 2.0, 3.0 };

    topography.set_height( 0, 0, 4 ); // Bottom row to 4
    topography.set_height( 1, 0, 4 ); // Bottom row to 4
    topography.set_height( 0, 1, 1 ); // Top row to 1
    topography.set_height( 1, 1, 1 ); // Top row to 1

    Flowy::Vector2 coord = { 1.5, 3 };
    // Expected values
    auto height_expected = 0.5 * ( topography.height_data( 0, 0 ) + topography.height_data( 0, 1 ) );
    Flowy::Vector2 slope_expected
        = { 0.0, -( topography.height_data( 0, 1 ) - topography.height_data( 0, 0 ) ) / topography.cell_size() };

    auto [height, slope] = topography.height_and_slope( coord );

    INFO( fmt::format( " height = {}\n", height ) );
    INFO( fmt::format( " slope = {}\n", fmt::streamed( slope ) ) );

    REQUIRE( height == height_expected );
    REQUIRE( slope == slope_expected );
}
