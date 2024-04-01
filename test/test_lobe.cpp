#include "catch2/matchers/catch_matchers.hpp"
#include "definitions.hpp"
#include "fmt/core.h"
#include "lobe.hpp"
#include "math.hpp"
#include "simulation.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>
#include <vector>

TEST_CASE( "Testing if a given point is inside a lobe or not", "[is_point_in_lobe]" )
{
    auto my_lobe            = Flowtastic::Lobe();
    my_lobe.azimuthal_angle = 0;        // The azimuthal angle of the major semi-axis with respect to the x-axis
    my_lobe.center          = { 0, 0 }; // The center of the ellipse
    my_lobe.semi_axes       = { 3, 1 };

    Flowtastic::Vector2 test_point_in  = { 2.0, 0.5 };
    Flowtastic::Vector2 test_point_out = { 1.0, 4.0 };

    REQUIRE( my_lobe.is_point_in_lobe( test_point_in ) );
    REQUIRE( !my_lobe.is_point_in_lobe( test_point_out ) );
}

TEST_CASE( "Test the bounding box of the lobe", "[lobe_bbox]" )
{
    auto my_lobe            = Flowtastic::Lobe();
    my_lobe.azimuthal_angle = Flowtastic::Math::pi / 2.0; // The azimuthal angle; major axis is along y-axis
    my_lobe.center          = { 0.5, 0 };                 // The center of the ellipse
    my_lobe.semi_axes       = { 2, 1 };
    // Expected bbox coordinates
    auto bbox_coord_expected
        = std::vector<Flowtastic::Vector2>{ { 1.5, 2.0 }, { -0.5, 2.0 }, { -0.5, -2.0 }, { 1.5, -2.0 } };

    auto bbox_coord = my_lobe.bounding_box();

    for( size_t i = 0; i < bbox_coord.size(); i++ )
    {
        for( size_t j = 0; j < 2; j++ )
        {
            REQUIRE_THAT( bbox_coord[i][j], Catch::Matchers::WithinRel( bbox_coord_expected[i][j] ) );
        }
    }
}