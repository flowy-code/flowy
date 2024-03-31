#include "catch2/matchers/catch_matchers.hpp"
#include "definitions.hpp"
#include "topography.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <xtensor/xio.hpp>

TEST_CASE( "height_and_slope_test", "[topography]" )
{
    auto topography = Flowtastic::Topography();
    // x and y axes have the usual meaning: we assume that the y axis is already "flipped" from the ASC file default
    // (top to down)
    topography.height_data = { { 4.0, 1.0 }, { 4.0, 1.0 } };
    topography.x_data      = { 1.0, 2.0 };
    topography.y_data      = { 2.0, 3.0 };

    Flowtastic::Vector2 coord = { 1.5, 2.5 };
    // Expected values
    auto height_expected = 0.5 * ( topography.height_data( 0, 0 ) + topography.height_data( 0, 1 ) );
    Flowtastic::Vector2 slope_expected
        = { ( topography.height_data( 0, 1 ) - topography.height_data( 0, 0 ) ) / topography.cell_size(), 0.0 };

    auto [height, slope] = topography.height_and_slope( coord );

    REQUIRE( height == height_expected );
    REQUIRE( slope == slope_expected );

    topography.height_data = { { 2.0, 1.0 }, { 3.0, 2.0 } };
    topography.x_data      = { 1.0, 2.0 };
    topography.y_data      = { 2.0, 3.0 };

    slope_expected            = { -1.0, 1.0 };
    std::tie( height, slope ) = topography.height_and_slope( coord );
    REQUIRE( slope == slope_expected );
}