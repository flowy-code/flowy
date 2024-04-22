// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/asc_file.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/topography_file.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>
#include <xtensor/xio.hpp>

TEST_CASE( "asc_file_test", "[asc]" )
{
    namespace fs = std::filesystem;

    auto proj_root_path = fs::current_path();
    auto asc_file_path  = proj_root_path / fs::path( "test/res/asc/file.asc" );

    auto asc_file = Flowy::AscFile( asc_file_path );

    INFO( fmt::format( "data = {}\n", fmt::streamed( asc_file.data ) ) );

    // NOTE: that the order of rows in the asc file is opposite to the order of rows in the data
    Flowy::MatrixX data_expected = { { 2.0, 0.5512 }, { 1.2, -9999 }, { 3, 2 } };

    REQUIRE( asc_file.data == data_expected );

    double xllcorner_expected     = 1e+02;
    double yllcorner_expected     = 2e+03;
    int cellsize_expected         = 100;
    double no_data_value_expected = -9999;

    Flowy::VectorX x_data_expected
        = { xllcorner_expected, xllcorner_expected + cellsize_expected, xllcorner_expected + 2.0 * cellsize_expected };
    Flowy::VectorX y_data_expected = { yllcorner_expected, yllcorner_expected + cellsize_expected };

    REQUIRE( xllcorner_expected == asc_file.lower_left_corner()[0] );
    REQUIRE( yllcorner_expected == asc_file.lower_left_corner()[1] );
    REQUIRE( cellsize_expected == asc_file.cell_size() );
    REQUIRE( no_data_value_expected == asc_file.no_data_value );
    REQUIRE( x_data_expected == asc_file.x_data );
    REQUIRE( y_data_expected == asc_file.y_data );
}

TEST_CASE( "asc_file_test_with_crop", "[asc_crop]" )
{
    namespace fs = std::filesystem;

    auto proj_root_path = fs::current_path();
    auto asc_file_path  = proj_root_path / fs::path( "test/res/asc/file_crop.asc" );

    Flowy::TopographyCrop crop{};

    crop.x_min = 1.2;
    crop.x_max = 8.2;
    crop.y_min = 1.2;
    crop.y_max = 5.2;

    auto asc_file = Flowy::AscFile( asc_file_path, crop );

    INFO( fmt::format( "data = {}\n", fmt::streamed( asc_file.data ) ) );
    INFO( fmt::format( "x_data = {}\n", fmt::streamed( asc_file.x_data ) ) );
    INFO( fmt::format( "y_data = {}\n", fmt::streamed( asc_file.y_data ) ) );
}
