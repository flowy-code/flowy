#include "definitions.hpp"
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

TEST_CASE( "Testing if a given point is inside a lobe or not", "[is_point_in_lobe]" )
{
    auto lobe            = Flowtastic::Lobe();
    lobe.azimuthal_angle = 0;        // The azimuthal angle of the major semi-axis with respect to the x-axis
    lobe.center          = { 0, 0 }; // The center of the ellipse
    lobe.semi_axes       = { 3, 1 };

    Flowtastic::Vector2 test_point_in  = { 2.0, 0.5 };
    Flowtastic::Vector2 test_point_out = { 1.0, 4.0 };

    REQUIRE( lobe.is_point_in_lobe( test_point_in ) );
    REQUIRE( !lobe.is_point_in_lobe( test_point_out ) );
}