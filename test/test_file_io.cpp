// GPL v3 License
// Copyright 2023--present Flowy developers
#include <netcdf.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include "flowy/include/asc_file.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/netcdf_file.hpp"
#include "flowy/include/topography.hpp"
#include "flowy/include/topography_file.hpp"
#include "fmt/ostream.h"
#include "xtensor/xbuilder.hpp"
#include <filesystem>

void test_read_write(
    const Flowy::Topography & topography, auto & write_cb, const auto & read_cb, auto & check_cb,
    const std::filesystem::path & path )
{
    INFO( path.string() );
    // Output path expects no extension, so we remove it
    auto path_write = path;
    path_write.replace_extension();
    write_cb( topography, path_write );
    auto file = read_cb( path );
    check_cb( topography, file );
}

TEST_CASE( "file_io", "[netcdf]" )
{
    namespace fs        = std::filesystem;
    auto proj_root_path = fs::current_path();

    auto output_folder = proj_root_path / "test" / "file_io";

    // Always start from a clean slate
    fs::remove_all( output_folder );
    fs::create_directories( output_folder );

    // Create a topography to test on
    Flowy::VectorX x_data      = xt::arange<double>( 0, 10, 1.0 );
    Flowy::VectorX y_data      = xt::arange<double>( 0, 2.0, 1.0 );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    for( size_t i = 0; i < x_data.size(); i++ )
    {
        for( size_t j = 0; j < y_data.size(); j++ )
        {
            height_data( i, j ) = x_data[i];
        }
    }

    double NO_DATA_VALUE = -1;

    auto topography   = Flowy::Topography( height_data, x_data, y_data, NO_DATA_VALUE );
    topography.hazard = xt::ones_like( height_data );

    auto write_asc_height = [&]( const Flowy::Topography & topography, const fs::path & path )
    {
        auto file = Flowy::AscFile( topography, Flowy::OutputQuantitiy::Height );
        file.save( path );
    };

    auto write_asc_hazard = [&]( const Flowy::Topography & topography, const fs::path & path )
    {
        auto file = Flowy::AscFile( topography, Flowy::OutputQuantitiy::Hazard );
        file.save( path );
    };

    auto read_asc            = [&]( const fs::path & path ) { return Flowy::AscFile( path ); };
    auto write_netcdf_height = [&]( const Flowy::Topography & topography, const fs::path & path )
    {
        auto file      = Flowy::NetCDFFile( topography, Flowy::OutputQuantitiy::Height );
        file.data_type = Flowy::StorageDataType::Double;
        file.save( path );
    };

    auto write_netcdf_height_float = [&]( const Flowy::Topography & topography, const fs::path & path )
    {
        auto file      = Flowy::NetCDFFile( topography, Flowy::OutputQuantitiy::Height );
        file.data_type = Flowy::StorageDataType::Float;
        file.save( path );
    };

    auto write_netcdf_height_short = [&]( const Flowy::Topography & topography, const fs::path & path )
    {
        auto file      = Flowy::NetCDFFile( topography, Flowy::OutputQuantitiy::Height );
        file.data_type = Flowy::StorageDataType::Short;
        file.save( path );
    };

    auto write_netcdf_hazard = [&]( const Flowy::Topography & topography, const fs::path & path )
    {
        auto file = Flowy::NetCDFFile( topography, Flowy::OutputQuantitiy::Hazard );
        file.save( path );
    };

    auto read_netcdf = [&]( const fs::path & path ) { return Flowy::NetCDFFile( path ); };

    auto check_height = [&]( const auto & topography, auto & file )
    {
        INFO( fmt::format( "file.x_data = {}\n", fmt::streamed( file.x_data ) ) );
        INFO( fmt::format( "file.y_data = {}\n", fmt::streamed( file.y_data ) ) );
        INFO( fmt::format( "top.x_data = {}\n", fmt::streamed( topography.x_data ) ) );
        INFO( fmt::format( "top.y_data = {}\n", fmt::streamed( topography.y_data ) ) );
        INFO( fmt::format( "file.data = {}\n", fmt::streamed( file.data ) ) );
        INFO( fmt::format( "top.height_data = {}\n", fmt::streamed( topography.height_data ) ) );
        REQUIRE( file.no_data_value == NO_DATA_VALUE );
        REQUIRE( xt::isclose( file.x_data, topography.x_data )() );
        REQUIRE( xt::isclose( file.y_data, topography.y_data )() );
        REQUIRE( xt::isclose( file.data, topography.height_data, 1e-4, 1e-4 )() );
    };

    auto check_hazard = [&]( const auto & topography, auto & file )
    {
        INFO( fmt::format( "file.x_data = {}\n", fmt::streamed( file.x_data ) ) );
        INFO( fmt::format( "file.y_data = {}\n", fmt::streamed( file.y_data ) ) );
        INFO( fmt::format( "top.x_data = {}\n", fmt::streamed( topography.x_data ) ) );
        INFO( fmt::format( "top.y_data = {}\n", fmt::streamed( topography.y_data ) ) );
        INFO( fmt::format( "file.data = {}\n", fmt::streamed( file.data ) ) );
        INFO( fmt::format( "top.hazard = {}\n", fmt::streamed( topography.hazard ) ) );
        REQUIRE( file.no_data_value == NO_DATA_VALUE );
        REQUIRE( xt::isclose( file.x_data, topography.x_data )() );
        REQUIRE( xt::isclose( file.y_data, topography.y_data )() );
        REQUIRE( xt::isclose( file.data, topography.hazard, 1e-4, 1e-4 )() );
    };

    test_read_write( topography, write_netcdf_height, read_netcdf, check_height, output_folder / "height.nc" );
    test_read_write(
        topography, write_netcdf_height_float, read_netcdf, check_height, output_folder / "height_float.nc" );
    test_read_write(
        topography, write_netcdf_height_short, read_netcdf, check_height, output_folder / "height_short.nc" );

    test_read_write( topography, write_asc_height, read_asc, check_height, output_folder / "height.asc" );
    test_read_write( topography, write_netcdf_hazard, read_netcdf, check_hazard, output_folder / "hazard.nc" );
    test_read_write( topography, write_asc_hazard, read_asc, check_hazard, output_folder / "hazard.asc" );
}

TEST_CASE( "file_test_with_crop", "[crop]" )
{
    namespace fs = std::filesystem;

    auto proj_root_path = fs::current_path();
    auto asc_file_path  = proj_root_path / fs::path( "test/res/asc/file_crop.asc" );

    Flowy::TopographyCrop crop{};

    crop.x_min = 2.2;
    crop.x_max = 8.2;
    crop.y_min = 1.2;
    crop.y_max = 5.2;

    // Test that the crop works when an .asc file is read in and used in the constructor of AscFile
    auto asc_file = Flowy::AscFile( asc_file_path, crop );

    INFO( fmt::format( "asc data = {}\n", fmt::streamed( asc_file.data ) ) );
    INFO( fmt::format( "asc x_data = {}\n", fmt::streamed( asc_file.x_data ) ) );
    INFO( fmt::format( "asc y_data = {}\n", fmt::streamed( asc_file.y_data ) ) );

    REQUIRE( asc_file.x_data[0] <= crop.x_min );
    REQUIRE( asc_file.x_data.periodic( -1 ) >= crop.x_max );
    REQUIRE( asc_file.y_data[0] <= crop.x_min );
    REQUIRE( asc_file.y_data.periodic( -1 ) >= crop.y_max );
    REQUIRE( asc_file.data.shape() == xt::shape( { asc_file.x_data.size(), asc_file.y_data.size() } ) );
}

TEST_CASE( "file_crop_to_content", "[crop_to_content]" )
{
    namespace fs = std::filesystem;

    auto proj_root_path   = fs::current_path();
    auto asc_file_path    = proj_root_path / fs::path( "test/res/asc/before_crop_to_content.asc" );
    auto output_file_path = proj_root_path / fs::path( "test/res/asc/after_crop_to_content" );

    // Read in the topography from the ASCII file, with no crop
    auto asc_file = Flowy::AscFile( asc_file_path );

    INFO( fmt::format( "Data from file before crop to content\n,data = {}\n", fmt::streamed( asc_file.data ) ) );
    INFO( fmt::format( "x_data = {}\n", fmt::streamed( asc_file.x_data ) ) );
    INFO( fmt::format( "y_data = {}\n", fmt::streamed( asc_file.y_data ) ) );

    double max_before = xt::amax( asc_file.data )();
    asc_file.crop_to_content();
    double max_after = xt::amax( asc_file.data )();

    REQUIRE( max_before == max_after );

    INFO( fmt::format( "Data from file after crop to content\n,data = {}\n", fmt::streamed( asc_file.data ) ) );
    INFO( fmt::format( "x_data = {}\n", fmt::streamed( asc_file.x_data ) ) );
    INFO( fmt::format( "y_data = {}\n", fmt::streamed( asc_file.y_data ) ) );

    asc_file.save( output_file_path );
}
