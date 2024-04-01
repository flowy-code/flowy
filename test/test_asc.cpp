#include "asc_file.hpp"
#include "definitions.hpp"
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

    auto asc_file = Flowtastic::AscFile( asc_file_path );

    fmt::print( "data = {}\n", fmt::streamed( asc_file.height_data ) );

    // NOTE: that the order of rows in the asc file is opposite to the order of rows in the height_data
    Flowtastic::MatrixX height_data_expected = { { 2.0, 1.2, 3 }, { 0.5512, -9999, 2 } };

    REQUIRE( asc_file.height_data == height_data_expected );

    double xllcorner_expected     = 1e+02;
    double yllcorner_expected     = 2e+03;
    int cellsize_expected         = 100;
    double no_data_value_expected = -9999;

    Flowtastic::VectorX x_data_expected
        = { xllcorner_expected, xllcorner_expected + cellsize_expected, xllcorner_expected + 2.0 * cellsize_expected };
    Flowtastic::VectorX y_data_expected = { yllcorner_expected, yllcorner_expected + cellsize_expected };

    REQUIRE( xllcorner_expected == asc_file.lower_left_corner[0] );
    REQUIRE( yllcorner_expected == asc_file.lower_left_corner[1] );
    REQUIRE( cellsize_expected == asc_file.cell_size );
    REQUIRE( no_data_value_expected == asc_file.no_data_value );
    REQUIRE( x_data_expected == asc_file.x_data );
    REQUIRE( y_data_expected == asc_file.y_data );
}