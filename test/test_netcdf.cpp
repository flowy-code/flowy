// GPL v3 License
// Copyright 2023--present Flowy developers
#include <netcdf.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include "flowy/include/definitions.hpp"

#define NC_EXAMPLE_ERROR 2 /* This is the exit code for failure. */

TEST_CASE( "netcdf_test", "[netcdf]" )
{
    Flowy::VectorX x_data      = xt::linspace<double>( 0, 10, 10, true );
    Flowy::VectorX y_data      = xt::linspace<double>( 0, 1.0, 2, true );
    Flowy::MatrixX height_data = xt::zeros<double>( { x_data.size(), y_data.size() } );

    for( size_t i = 0; i < x_data.size(); i++ )
    {
        for( size_t j = 0; j < y_data.size(); j++ )
        {
            height_data( i, j ) = x_data[i];
        }
    }

    // height_data = xt::transpose( xt::flip( height_data, 0 ) )

    height_data( 0, 0 ) = -9999;

    double nodata_value = -9999;

    std::string filename = "test.nc";
    std::string unit     = "meters";

    std::string missing_value = "missing_value";
    std::string fill_value    = "_FillValue";

    int ncid;
    int dimids[2];
    int numdims = height_data.shape().size();

    nc_create( filename.c_str(), NC_CLOBBER, &ncid );

    nc_def_dim( ncid, "y", height_data.shape()[0], dimids );
    nc_def_dim( ncid, "x", height_data.shape()[1], &dimids[1] );

    // header for x variable
    int varid_x;
    int dimid_x;
    nc_def_var( ncid, "x", NC_DOUBLE, 1, &dimid_x, &varid_x );
    nc_put_att_text( ncid, varid_x, "units", unit.size(), unit.c_str() );

    // header for y variable
    int varid_y;
    int dimid_y;
    nc_def_var( ncid, "y", NC_DOUBLE, 1, &dimid_y, &varid_y );
    nc_put_att_text( ncid, varid_y, "units", unit.size(), unit.c_str() );

    // header for elevation variable
    int varid_elevation;
    nc_def_var( ncid, "elevation", NC_DOUBLE, numdims, dimids, &varid_elevation );
    nc_put_att_text( ncid, varid_elevation, "units", unit.size(), unit.c_str() );

    nc_def_var_fill( ncid, varid_elevation, 0, &nodata_value );

    nc_enddef( ncid );

    // Write elevation data
    nc_put_var_double( ncid, varid_elevation, (double *)height_data.data() );
    // Write x data
    nc_put_var_double( ncid, varid_x, (double *)x_data.data() );
    // Write y data
    nc_put_var_double( ncid, varid_y, (double *)y_data.data() );

    nc_close( ncid );
}

TEST_CASE( "asc_file_test_with_crop", "[asc_crop]" )
{
    // namespace fs = std::filesystem;

    // auto proj_root_path = fs::current_path();
    // auto asc_file_path  = proj_root_path / fs::path( "test/res/asc/file_crop.asc" );

    // Flowy::AscCrop crop{};

    // crop.x_min = 1.2;
    // crop.x_max = 8.2;
    // crop.y_min = 1.2;
    // crop.y_max = 5.2;

    // auto asc_file = Flowy::AscFile( asc_file_path, crop );

    // INFO( fmt::format( "data = {}\n", fmt::streamed( asc_file.height_data ) ) );
    // INFO( fmt::format( "x_data = {}\n", fmt::streamed( asc_file.x_data ) ) );
    // INFO( fmt::format( "y_data = {}\n", fmt::streamed( asc_file.y_data ) ) );
}
