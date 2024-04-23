// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/netcdf_file.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/dump_csv.hpp"
#include "xtensor/xmanipulation.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xview.hpp"
#include <fmt/format.h>
#include <netcdf.h>
#include <algorithm>
#include <fstream>
#include <limits>
#include <stdexcept>

namespace Flowy
{

void NetCDFFile::determine_scale_and_offset()
{
    double data_max = xt::amax( data )();
    double data_min = xt::amin( data )();

    auto int_range = static_cast<double>( std::numeric_limits<short>::max() )
                     - static_cast<double>( std::numeric_limits<short>::min() );

    scale_factor = ( data_max - data_min ) / int_range;
    add_offset   = data_min - std::numeric_limits<short>::min() * scale_factor.value();
}

void CHECK( int retval )
{
    if( static_cast<bool>( retval ) )
    {
        throw std::runtime_error( fmt::format( "NetCDF error {}\n", nc_strerror( retval ) ) );
    }
}

// Wrap the NetCDF file handle in RAII struct
struct FileHandle
{
    int ncid{};
    FileHandle( const std::filesystem::path & path )
    {
        CHECK( nc_create( path.string().c_str(), NC_NETCDF4, &ncid ) );
    }

    ~FileHandle()
    {
        nc_close( ncid );
    }
};

void NetCDFFile::save( const std::filesystem::path & path )
{

    // define data type according to desired output type
    int dtype{};
    if( data_type == StorageDataType::Double )
    {
        dtype = NC_DOUBLE;
    }
    else if( data_type == StorageDataType::Float )
    {
        dtype = NC_FLOAT;
    }
    else
    {
        determine_scale_and_offset();
        no_data_value = std::round( ( no_data_value - add_offset.value() ) / scale_factor.value() );
        dtype         = NC_SHORT;
    }

    const std::string unit = "meters";

    std::string missing_value = "missing_value";
    std::string fill_value    = "_FillValue";
    MatrixX data_out          = xt::flip( xt::transpose( data ), 0 );

    int dimids[2];
    int numdims = data_out.shape().size();

    auto file_handle = FileHandle( path );
    int ncid         = file_handle.ncid;

    CHECK( nc_def_dim( ncid, "x", data_out.shape()[0], dimids ) );
    CHECK( nc_def_dim( ncid, "y", data_out.shape()[1], &dimids[1] ) );

    // header for y variable (we save x_data here)
    int varid_y{};
    int dimid_y{};
    CHECK( nc_def_var( ncid, "y", NC_DOUBLE, 1, &dimid_y, &varid_y ) );
    CHECK( nc_put_att_text( ncid, varid_y, "units", unit.size(), unit.c_str() ) );

    // header for x variable (we save y_data here)
    int varid_x{};
    int dimid_x{};
    CHECK( nc_def_var( ncid, "x", NC_DOUBLE, 1, &dimid_x, &varid_x ) );
    CHECK( nc_put_att_text( ncid, varid_x, "units", unit.size(), unit.c_str() ) );

    // header for elevation variable
    int varid_elevation{};

    CHECK( nc_def_var( ncid, "elevation", dtype, numdims, dimids, &varid_elevation ) );
    CHECK( nc_put_att_text( ncid, varid_elevation, "units", unit.size(), unit.c_str() ) );

    CHECK( nc_def_var_fill( ncid, varid_elevation, 0, &no_data_value ) );

    if( compression )
    {
        CHECK( nc_def_var_deflate( ncid, varid_elevation, static_cast<int>( shuffle ), 1, compression_level ) );
    }

    if( scale_factor.has_value() )
    {
        double temp_scale_factor = scale_factor.value();
        CHECK( nc_put_att_double( ncid, varid_elevation, "scale_factor", NC_DOUBLE, 1, &temp_scale_factor ) );
    }

    if( add_offset.has_value() )
    {
        double temp_add_offset = add_offset.value();
        CHECK( nc_put_att_double( ncid, varid_elevation, "add_offset", NC_DOUBLE, 1, &temp_add_offset ) );
    }

    // end of header definition
    CHECK( nc_enddef( ncid ) );

    // Write elevation data
    if( data_type == StorageDataType::Double )
    {
        CHECK( nc_put_var_double( ncid, varid_elevation, static_cast<double *>( data_out.data() ) ) );
    }
    else if( data_type == StorageDataType::Float )
    {
        xt::xtensor<float, 2> data_ = data_out;
        CHECK( nc_put_var_float( ncid, varid_elevation, static_cast<float *>( data_.data() ) ) );
    }
    else if( data_type == StorageDataType::Short )
    {
        xt::xtensor<short, 2> data_ = xt::round( ( data_out - add_offset.value() ) / scale_factor.value() );
        CHECK( nc_put_var_short( ncid, varid_elevation, static_cast<short *>( data_.data() ) ) );
    }

    // Write x data
    // Since we transpose the height data, we have to write our x_data to the netcdf "y" attribute
    // QGIS also expects the coordinates of pixel centers, so we shift by 0.5 * cell_size
    VectorX x_data_out = x_data + 0.5 * cell_size();
    CHECK( nc_put_var_double( ncid, varid_y, static_cast<double *>( x_data_out.data() ) ) );

    // Write y data
    // Since we transpose the height data and flip the y values, we have to write our y_data to the netcdf "x"
    // attribute and flip it QGIS also expects the coordinates of pixel centers, so we shift by 0.5 * cell_size
    VectorX y_data_out = xt::flip( y_data ) + 0.5 * cell_size();
    CHECK( nc_put_var_double( ncid, varid_x, static_cast<double *>( y_data_out.data() ) ) );
}

} // namespace Flowy
