// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/math.hpp"
#include "xtensor/xmath.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <optional>

TEST_CASE( "Testing if a given point is inside a lobe or not", "[is_point_in_lobe]" )
{
    auto my_lobe = Flowy::Lobe();
    my_lobe.set_azimuthal_angle( 0 ); // The azimuthal angle of the major semi-axis with respect to the x-axis
    my_lobe.center    = { 0, 0 };     // The center of the ellipse
    my_lobe.semi_axes = { 3, 1 };

    Flowy::Vector2 test_point_in  = { 2.0, 0.5 };
    Flowy::Vector2 test_point_out = { 1.0, 4.0 };

    REQUIRE( my_lobe.is_point_in_lobe( test_point_in ) );
    REQUIRE( !my_lobe.is_point_in_lobe( test_point_out ) );
}

void check_intersection(
    const Flowy::Lobe & lobe, const Flowy::Vector2 & x1, const Flowy::Vector2 & x2,
    const std::optional<std::array<Flowy::Vector2, 2>> & expected )
{
    auto intersection = lobe.line_segment_intersects( x1, x2 );

    INFO( fmt::format( "x1 = {}\n", fmt::streamed( x1 ) ) );
    INFO( fmt::format( "x2 = {}\n", fmt::streamed( x2 ) ) );
    INFO( fmt::format( "intersection = {}\n", intersection ) );
    INFO( fmt::format( "expected = {}\n", expected ) );

    if( intersection.has_value() )
    {
        REQUIRE( xt::isclose( expected.value()[0], expected.value()[0] )() );
        REQUIRE( xt::isclose( expected.value()[1], expected.value()[1] )() );
    }
    else
    {
        REQUIRE( !expected.has_value() );
    }
}

TEST_CASE( "Testing if a line segment intersects a lobe", "[line_segment_intersects]" )
{
    using Flowy::Vector2;

    auto my_lobe = Flowy::Lobe();
    my_lobe.set_azimuthal_angle( 0 ); // The azimuthal angle of the major semi-axis with respect to the x-axis
    my_lobe.center    = { 1, 2 };     // The center of the ellipse
    my_lobe.semi_axes = { 2, 1 };
    my_lobe.set_azimuthal_angle( Flowy::Math::pi / 4 );

    // These values have been tested in flowpy
    Vector2 x1                      = { 1.0, 2.0 };
    Vector2 x2                      = { 2, 2 };
    Vector2 p1                      = { 1.0, 2.0 };
    Vector2 p2                      = { 2.0, 2.0 };
    std::array<Vector2, 2> expected = { p1, p2 };
    check_intersection( my_lobe, x1, x2, expected );

    x1       = { 1.0, 2.0 };
    x2       = { 4, 4 };
    p1       = { 1.0, 2.0 };
    p2       = { 2.575677194316671, 3.0504514628777804 };
    expected = { p1, p2 };
    check_intersection( my_lobe, x1, x2, expected );

    x1       = { -1, -1 };
    x2       = { 4, 4 };
    p1       = { 0.4999999999999999, 0.49999999999999956 };
    p2       = { 2.5000000000000013, 2.500000000000001 };
    expected = { p1, p2 };
    check_intersection( my_lobe, x1, x2, expected );

    x1 = { -1, -1 };
    x2 = { 4, -4 };
    check_intersection( my_lobe, x1, x2, std::nullopt );

    x1       = { 0, 2.5 };
    x2       = { 0.5, 2.5 };
    p1       = { 0.10000000000000009, 2.5 };
    p2       = { 0.5, 2.5 };
    expected = { p1, p2 };
    check_intersection( my_lobe, x1, x2, expected );

    x1 = { 0, 3.0 };
    x2 = { -0.5, 3.0 };
    check_intersection( my_lobe, x1, x2, std::nullopt );
}
