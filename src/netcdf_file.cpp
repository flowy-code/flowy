// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/netcdf_file.hpp"
#include "flowy/include/definitions.hpp"
#include "xtensor/xmanipulation.hpp"
#include "xtensor/xmath.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <netcdf.h>
#include <limits>
#include <stdexcept>

namespace Flowy
{

void NetCDFFile::determine_scale_and_offset()
{
    double data_max = xt::amax( data )();
    data_max        = std::max( data_max, no_data_value );

    double data_min = xt::amin( data )();
    data_min        = std::min( data_min, no_data_value );

    auto int_range = static_cast<double>( std::numeric_limits<int16_t>::max() )
                     - static_cast<double>( std::numeric_limits<int16_t>::min() );

    scale_factor = ( data_max - data_min ) / int_range;
    add_offset   = data_min - std::numeric_limits<int16_t>::min() * scale_factor.value();
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
    explicit FileHandle( const std::filesystem::path & path, bool read = false )
    {
        if( read )
        {
            CHECK( nc_open( path.string().c_str(), NC_NETCDF4, &ncid ) );
        }
        else
        {
            CHECK( nc_create( path.string().c_str(), NC_NETCDF4, &ncid ) );
        }
    }

    ~FileHandle()
    {
        nc_close( ncid );
    }
};

enum class XY
{
    X,
    Y
};

VectorX get_xy_from_netcdf( int ncid, XY xy )
{
    // get the varid
    int varid{};
    std::string name{};
    if( xy == XY::X )
    {
        name = "x";
        CHECK( nc_inq_varid( ncid, "x", &varid ) );
    }
    else
    {
        name = "y";
        CHECK( nc_inq_varid( ncid, "y", &varid ) );
    }

    // read the attributes
    nc_type vartype{};
    int num_dims{};
    int dimids[NC_MAX_DIMS];
    int num_attrs{};

    CHECK( nc_inq_var( ncid, varid, nullptr, &vartype, &num_dims, dimids, &num_attrs ) );

    // compute the length
    size_t dim_lengths[NC_MAX_DIMS];

    for( int i = 0; i < num_dims; i++ )
    {
        CHECK( nc_inq_dimid( ncid, name.c_str(), &dimids[i] ) );
        CHECK( nc_inq_dimlen( ncid, dimids[i], &dim_lengths[i] ) );
    }

    size_t total_size = 1;
    for( int i = 0; i < num_dims; i++ )
    {
        total_size *= dim_lengths[i];
    }

    VectorX data = xt::zeros<double>( { total_size } );

    nc_get_var_double( ncid, varid, data.data() );

    return data;
}

MatrixX NetCDFFile::get_elevation_from_netcdf( int ncid )
{
    // get the varid
    int varid{};

    CHECK( nc_inq_varid( ncid, "elevation", &varid ) );

    // read the attributes
    nc_type vartype{};
    int num_dims{};
    int dimids[NC_MAX_DIMS];
    int num_attrs{};
    CHECK( nc_inq_var( ncid, varid, nullptr, &vartype, &num_dims, dimids, &num_attrs ) );

    // compute the length
    size_t dim_lengths[NC_MAX_DIMS];

    CHECK( nc_inq_dimid( ncid, "x", &dimids[0] ) );
    CHECK( nc_inq_dimid( ncid, "y", &dimids[1] ) );

    for( int i = 0; i < num_dims; i++ )
    {
        CHECK( nc_inq_dimlen( ncid, dimids[i], &dim_lengths[i] ) );
    }

    // Set the no_data value to the fill value
    CHECK( nc_get_att_double( ncid, varid, "_FillValue", &no_data_value ) );

    // This is the data we will return
    MatrixX data = xt::zeros<double>( { dim_lengths[0], dim_lengths[1] } );

    if( vartype == NC_FLOAT )
    {
        xt::xtensor<float, 2> data_read = xt::zeros<float>( { dim_lengths[0], dim_lengths[1] } );
        nc_get_var_float( ncid, varid, data_read.data() );
        data = data_read;
    }
    else if( vartype == NC_DOUBLE )
    {
        nc_get_var_double( ncid, varid, data.data() );
    }
    else if( vartype == NC_SHORT )
    {
        // Read the scale factor and add_offset attributes
        double scale_factor{};
        CHECK( nc_get_att_double( ncid, varid, "scale_factor", &scale_factor ) );

        double add_offset{};
        CHECK( nc_get_att_double( ncid, varid, "add_offset", &add_offset ) );

        xt::xtensor<int16_t, 2> data_read = xt::zeros<int16_t>( { dim_lengths[0], dim_lengths[1] } );
        nc_get_var_short( ncid, varid, data_read.data() );
        data = data_read;
        data *= scale_factor;
        data += add_offset;

        no_data_value *= scale_factor;
        no_data_value += add_offset;
    }

    // We need to undo the transformation we apply before saving
    data = xt::eval( xt::transpose( xt::flip( data, 0 ) ) );
    return data;
}

NetCDFFile::NetCDFFile( const std::filesystem::path & path, const std::optional<TopographyCrop> & crop )
{
    auto file_handle = FileHandle( path, true );
    int ncid         = file_handle.ncid;

    // Remember: Our x,y and their xy are flipped. Also theirs refers to the center of pixels, while ours refers to the
    // lower left corner
    y_data = get_xy_from_netcdf( ncid, XY::X );
    y_data = xt::flip( y_data );
    y_data -= 0.5 * cell_size();
    x_data = get_xy_from_netcdf( ncid, XY::Y );
    x_data -= 0.5 * cell_size();
    data = get_elevation_from_netcdf( ncid );

    // apply crop
    if( crop.has_value() )
        crop_topography( crop.value() );
}

void NetCDFFile::save( const std::filesystem::path & path_ )
{
    auto path = handle_suffix( path_ );

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
    CHECK( nc_def_var( ncid, "y", NC_DOUBLE, 1, &dimids[1], &varid_y ) );
    CHECK( nc_put_att_text( ncid, varid_y, "units", unit.size(), unit.c_str() ) );
    std::string att = "projection_y_coordinate";
    CHECK( nc_put_att_text( ncid, varid_y, "standard_name", att.size(), att.c_str() ) );

    // header for x variable (we save y_data here)
    int varid_x{};
    CHECK( nc_def_var( ncid, "x", NC_DOUBLE, 1, &dimids[0], &varid_x ) );
    CHECK( nc_put_att_text( ncid, varid_x, "units", unit.size(), unit.c_str() ) );
    att = "projection_x_coordinate";
    CHECK( nc_put_att_text( ncid, varid_x, "standard_name", att.size(), att.c_str() ) );

    // header for elevation variable
    int varid_elevation{};
    CHECK( nc_def_var( ncid, "elevation", dtype, numdims, dimids, &varid_elevation ) );
    CHECK( nc_put_att_text( ncid, varid_elevation, "units", unit.size(), unit.c_str() ) );

    // Set the fill value
    if( data_type == StorageDataType::Double )
    {
        CHECK( nc_def_var_fill( ncid, varid_elevation, 0, &no_data_value ) );
    }
    else if( data_type == StorageDataType::Float )
    {
        float no_data_value_temp = no_data_value;
        CHECK( nc_def_var_fill( ncid, varid_elevation, 0, &no_data_value_temp ) );
    }
    else if( data_type == StorageDataType::Short )
    {
        int16_t no_data_value_temp = no_data_value;
        CHECK( nc_def_var_fill( ncid, varid_elevation, 0, &no_data_value_temp ) );
    }

    // Write elevation data
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
        xt::xtensor<int16_t, 2> data_ = xt::round( ( data_out - add_offset.value() ) / scale_factor.value() );
        CHECK( nc_put_var_short( ncid, varid_elevation, static_cast<int16_t *>( data_.data() ) ) );
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
