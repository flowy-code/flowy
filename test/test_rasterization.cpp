// GPL v3 License
// Copyright 2023--present Flowy developers
#include "catch2/matchers/catch_matchers.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/math.hpp"
#include "flowy/include/topography.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <xtensor/xio.hpp>

TEST_CASE( "test_rasterization", "[topography]" )
{
    Flowy::VectorX x_data      = xt::linspace( 0.0, 40.0, 40 );
    Flowy::VectorX y_data      = xt::linspace( 0.0, 20.0, 20 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto lobe      = Flowy::Lobe();
    lobe.semi_axes = { 8.0, 2.0 };
    lobe.thickness = 20.0;
    lobe.set_azimuthal_angle( Flowy::Math::pi / 4 );
    lobe.center = { 20.0, 10.0 };

    auto topography = Flowy::Topography( height_data, x_data, y_data );
    auto lobe_cells = topography.get_cells_intersecting_lobe( lobe );

    std::vector<std::array<int, 2>> cells_intersected_expected = {
        { 13, 4 },  { 14, 4 },  { 15, 4 },  { 16, 4 },  { 13, 5 },  { 14, 5 },  { 16, 5 },  { 17, 5 },  { 18, 5 },
        { 14, 6 },  { 18, 6 },  { 19, 6 },  { 14, 7 },  { 15, 7 },  { 19, 7 },  { 20, 7 },  { 15, 8 },  { 16, 8 },
        { 20, 8 },  { 21, 8 },  { 16, 9 },  { 17, 9 },  { 21, 9 },  { 22, 9 },  { 17, 10 }, { 18, 10 }, { 22, 10 },
        { 23, 10 }, { 18, 11 }, { 19, 11 }, { 23, 11 }, { 24, 11 }, { 19, 12 }, { 20, 12 }, { 24, 12 }, { 20, 13 },
        { 21, 13 }, { 22, 13 }, { 24, 13 }, { 25, 13 }, { 22, 14 }, { 23, 14 }, { 24, 14 }, { 25, 14 }, { 24, 15 }
    };

    std::vector<std::array<int, 2>> cells_enclosed_expected
        = { { 15, 5 },  { 15, 6 },  { 16, 6 },  { 17, 6 },  { 16, 7 },  { 17, 7 },  { 18, 7 },  { 17, 8 },
            { 18, 8 },  { 19, 8 },  { 18, 9 },  { 19, 9 },  { 20, 9 },  { 19, 10 }, { 20, 10 }, { 21, 10 },
            { 20, 11 }, { 21, 11 }, { 22, 11 }, { 21, 12 }, { 22, 12 }, { 23, 12 }, { 23, 13 } };

    REQUIRE_THAT( cells_enclosed_expected, Catch::Matchers::UnorderedRangeEquals( lobe_cells.cells_enclosed ) );
    REQUIRE_THAT( cells_intersected_expected, Catch::Matchers::UnorderedRangeEquals( lobe_cells.cells_intersecting ) );
}

TEST_CASE( "test_add_lobe", "[topography]" )
{
    Flowy::VectorX x_data      = xt::linspace( 0.0, 40.0, 40 );
    Flowy::VectorX y_data      = xt::linspace( 0.0, 20.0, 20 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto lobe      = Flowy::Lobe();
    lobe.semi_axes = { 8.0, 2.0 };
    lobe.thickness = 20.0;
    lobe.set_azimuthal_angle( Flowy::Math::pi / 4 );
    lobe.center = { 20.0, 10.0 };

    auto topography = Flowy::Topography( height_data, x_data, y_data );
    auto lobe_cells = topography.get_cells_intersecting_lobe( lobe );

    auto volume_before = topography.volume();
    auto area_before   = topography.area( 0.0 );

    topography.add_lobe( lobe );

    auto volume_after = topography.volume();
    auto area_after   = topography.area( 0.0 );

    INFO( fmt::format(
        "volume_before = {}, volume_after = {}, diff = {}\n", volume_before, volume_after,
        volume_after - volume_before ) );
    INFO( fmt::format(
        "area_before = {}, area_after = {}, diff = {}\n", area_before, area_after, area_after - area_before ) );
    INFO( fmt::format( "True lobe area = {}\n", lobe.area() ) );
    INFO( fmt::format( "True lobe volume = {}\n", lobe.volume() ) );

    REQUIRE_THAT( volume_after - volume_before, Catch::Matchers::WithinRel( lobe.volume(), 0.05 ) );
}