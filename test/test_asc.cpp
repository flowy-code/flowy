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

    Flowtastic::MatrixX height_data_expected = { { 0.5512, -9999, 2 }, { 2.0, 1.2, 3 } };

    REQUIRE( asc_file.height_data == height_data_expected );

    double xllcorner_expected     = 2.701332e+05;
    double yllcorner_expected     = 2.123588e+06;
    int cellsize_expected         = 20;
    double no_data_value_expected = -9999;

    REQUIRE( xllcorner_expected == asc_file.lower_left_corner[0] );
    REQUIRE( yllcorner_expected == asc_file.lower_left_corner[1] );
    REQUIRE( cellsize_expected == asc_file.cell_size );
    REQUIRE( no_data_value_expected == asc_file.no_data_value );
}