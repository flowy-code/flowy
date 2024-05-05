// GPL v3 License
// Copyright 2023--present Flowy developers
#include "catch2/matchers/internal/catch_matchers_impl.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/math.hpp"
#include "flowy/include/topography.hpp"
#include "fmt/core.h"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xio.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/ranges.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <optional>
#include <vector>

TEST_CASE( "bounding_box", "[bounding_box]" )
{
    Flowy::VectorX x_data      = xt::arange<double>( 0.5, 19.5, 1.0 );
    Flowy::VectorX y_data      = xt::arange<double>( 4.5, 10.5, 1.0 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );

    Flowy::Vector2 point = { 11.4, 7.6 };

    int idx_point_x_expected  = 10;
    int idx_point_y_expected  = 3;
    double cell_size_expected = 1.0;

    auto point_indices = topography.locate_point( point );
    REQUIRE( idx_point_x_expected == point_indices[0] );
    REQUIRE( idx_point_y_expected == point_indices[1] );
    REQUIRE_THAT( topography.cell_size(), Catch::Matchers::WithinRel( cell_size_expected ) );

    INFO( fmt::format( "point.idx_x = {}, point.idx_y = {}\n", point_indices[0], point_indices[1] ) );

    // Note: this bounding box is chosen such that no clamping will occur
    auto bbox = topography.bounding_box( point, 2, 2 );
    INFO( fmt::format( "bbox.idx_x_higher = {}\n", bbox.idx_x_higher ) );
    INFO( fmt::format( "bbox.idx_x_lower = {}\n", bbox.idx_x_lower ) );
    INFO( fmt::format( "bbox.idx_y_lower = {}\n", bbox.idx_y_lower ) );
    INFO( fmt::format( "bbox.idx_y_higher = {}\n", bbox.idx_y_higher ) );

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
    Flowy::VectorX x_data      = xt::arange<double>( -2.0, 2.0, 1.0 );
    Flowy::VectorX y_data      = xt::arange<double>( -2.0, 2.0, 1.0 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } ); // not really needed here

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );

    Flowy::Lobe my_lobe;
    my_lobe.center    = { 0.0, 0.0 };
    my_lobe.semi_axes = { 1.95, 0.95 };
    my_lobe.set_azimuthal_angle( Flowy::Math::pi / 2.0 );

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

    auto lobe_cells = topography.get_cells_intersecting_lobe( my_lobe );

    REQUIRE_THAT( cell_indices_expected, Catch::Matchers::UnorderedRangeEquals( lobe_cells.cells_intersecting ) );
}

TEST_CASE( "test_compute_intersection", "[intersection]" )
{
    Flowy::VectorX x_data      = xt::arange<double>( -3, 3, 1.0 );
    Flowy::VectorX y_data      = xt::arange<double>( -3, 3, 1.0 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );

    Flowy::Lobe my_lobe;
    my_lobe.center    = { 0, 0 };
    my_lobe.semi_axes = { 1 - 1e-14, 1 - 1e-14 };
    my_lobe.set_azimuthal_angle( 0 );

    // clang-format off
    std::vector<std::array<int, 2>> cell_indices_expected = {
        { 2, 2 },
        { 2, 3 },
        { 3, 2 },
        { 3, 3 } };
    // clang-format on
    double expected_area_fraction = Flowy::Math::pi / 4.0;

    auto intersection_data = topography.compute_intersection( my_lobe, std::nullopt );

    std::vector<std::array<int, 2>> cell_indices{};
    for( size_t i = 0; i < cell_indices_expected.size(); i++ )
    {
        auto indices  = intersection_data[i].first;
        auto fraction = intersection_data[i].second;
        cell_indices.push_back( indices );
        INFO( fmt::format(
            " indices = {}, fraction = {}, fraction_expected {}\n", indices, fraction, expected_area_fraction ) );
        REQUIRE_THAT( fraction, Catch::Matchers::WithinRel( expected_area_fraction, 5e-2 ) );
    }

    REQUIRE_THAT( cell_indices_expected, Catch::Matchers::UnorderedRangeEquals( cell_indices ) );
}

TEST_CASE( "find_preliminary_budding_point", "[budding_point]" )
{
    Flowy::VectorX x_data      = xt::arange<double>( -2, 2, 1.0 );
    Flowy::VectorX y_data      = xt::arange<double>( -2, 2, 1.0 );
    Flowy::MatrixX height_data = 5.0 * xt::ones<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );

    topography.set_height( { 0.5, 0.5 }, -5.0 );

    Flowy::Lobe my_lobe;
    my_lobe.center    = { 0, 0 };
    my_lobe.semi_axes = { 0.8, 0.8 };
    my_lobe.set_azimuthal_angle( 0.0 );

    auto perimeter = my_lobe.rasterize_perimeter( 32 );
    for( auto & p : perimeter )
    {
        INFO( fmt::format( "p = {}\n", fmt::streamed( p ) ) );
        INFO( fmt::format( "height = {}\n\n", topography.height_and_slope( p ).first ) );
    }

    Flowy::Vector2 budding_point = topography.find_preliminary_budding_point( my_lobe, 32 );

    INFO( fmt::format( "budding_point  = {}", budding_point ) );

    // The budding point should be on the diagonal
    REQUIRE_THAT( budding_point[0], Catch::Matchers::WithinRel( budding_point[1] ) );
}

TEST_CASE( "test_finite_difference_slope", "[finite_difference_slope]" )
{
    Flowy::VectorX x_data      = xt::arange<double>( 0.0, 10.0, 1.0 );
    Flowy::VectorX y_data      = xt::arange<double>( 0.0, 1.1, 1.0 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    for( size_t i = 0; i < x_data.size(); i++ )
    {
        for( size_t j = 0; j < y_data.size(); j++ )
        {
            height_data( i, j ) = x_data[i];
        }
    }

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );

    const Flowy::Vector2 point1 = { 1.0, 0.5 };
    const Flowy::Vector2 point2 = { 3.0, 0.5 };
    double slope                = topography.slope_between_points( point1, point2 );
    double expected_slope       = 0.0;

    INFO( fmt::format( "slope={}\n", slope ) );
    // This would just be zero because of the default min_slope_drop=0.0
    REQUIRE_THAT( slope, Catch::Matchers::WithinRel( expected_slope ) );

    slope          = topography.slope_between_points( point2, point1 );
    expected_slope = 1.0;
    INFO( fmt::format( "slope={}\n", slope ) );
    REQUIRE_THAT( slope, Catch::Matchers::WithinRel( expected_slope ) );

    // Test that if you give it a minimum slope_drop it is set to that
    slope          = topography.slope_between_points( point1, point2, std::nullopt );
    expected_slope = -1.0;
    INFO( fmt::format( "slope={}\n", slope ) );
    REQUIRE_THAT( slope, Catch::Matchers::WithinRel( expected_slope ) );
}
TEST_CASE( "test_volume_correction", "[volume_correction]" )
{
    Flowy::VectorX x_data      = xt::arange<double>( -3, 3, 1.2 );
    Flowy::VectorX y_data      = xt::arange<double>( -3, 3, 1.2 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );

    Flowy::Lobe my_lobe;
    my_lobe.center    = { 0, 0 };
    my_lobe.semi_axes = { 2.0, 2.0 };
    my_lobe.thickness = 1.0;
    my_lobe.set_azimuthal_angle( 0 );

    auto intersection_data = topography.compute_intersection( my_lobe, std::nullopt );

    // Get the volume added from the topography
    auto volume_before_lobe = topography.volume();
    topography.add_lobe( my_lobe, true );
    auto volume_added = topography.volume() - volume_before_lobe;
    INFO( fmt::format( "Volume added = {}", volume_added ) );
    // True volume and true area of the lobe
    auto true_volume = my_lobe.volume();
    INFO( fmt::format( "True volume of the lobe = {}", true_volume ) );
    REQUIRE_THAT( volume_added, Catch::Matchers::WithinRel( true_volume ) );
}

TEST_CASE( "test_adding_topography", "[add_to_topography]" )
{
    Flowy::VectorX x_data      = xt::arange<double>( -3, 3, 1.2 );
    Flowy::VectorX y_data      = xt::arange<double>( -3, 3, 1.2 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    auto topography = Flowy::Topography( height_data, x_data, y_data, Flowy::DEFAULT_NO_DATA_VALUE_HEIGHT );
    topography.height_data( 0, 0 ) = 2;

    auto topography_to_add = topography;

    // Add one to the diagonal (this relies on x_data and y_data having the same size)
    topography_to_add.height_data += xt::diag( xt::ones<double>( { x_data.size() } ) );

    Flowy::MatrixX height_data_expected = topography.height_data + topography_to_add.height_data;

    topography.add_to_topography( topography_to_add );

    REQUIRE( xt::isclose( topography.height_data, height_data_expected )() );
}
