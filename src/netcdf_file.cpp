// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/netcdf_file.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/dump_csv.hpp"
#include "xtensor/xmanipulation.hpp"
#include "xtensor/xslice.hpp"
#include "xtensor/xview.hpp"
#include <fmt/format.h>
#include <netcdf.h>
#include <fstream>

namespace Flowy
{

// NetCDFFile::NetCDFFile( const std::filesystem::path & path, std::optional<AscCrop> crop )
// {
//
// }

void NetCDFFile::save( const std::filesystem::path & path )
{
    const std::string unit = "meters";

    std::string missing_value = "missing_value";
    std::string fill_value    = "_FillValue";

    MatrixX height_data_out = xt::flip( xt::transpose( height_data ), 0 );

    int ncid;
    int dimids[2];
    int numdims = height_data_out.shape().size();

    nc_create( path.string().c_str(), NC_CLOBBER, &ncid );

    nc_def_dim( ncid, "x", height_data_out.shape()[0], dimids );
    nc_def_dim( ncid, "y", height_data_out.shape()[1], &dimids[1] );

    // header for y variable (we save x_data here)
    int varid_y{};
    int dimid_y{};
    nc_def_var( ncid, "y", NC_DOUBLE, 1, &dimid_y, &varid_y );
    nc_put_att_text( ncid, varid_y, "units", unit.size(), unit.c_str() );

    // header for x variable (we save y_data here)
    int varid_x{};
    int dimid_x{};
    nc_def_var( ncid, "x", NC_DOUBLE, 1, &dimid_x, &varid_x );
    nc_put_att_text( ncid, varid_x, "units", unit.size(), unit.c_str() );

    // header for elevation variable
    int varid_elevation{};
    nc_def_var( ncid, "elevation", NC_DOUBLE, numdims, dimids, &varid_elevation );
    nc_put_att_text( ncid, varid_elevation, "units", unit.size(), unit.c_str() );

    // Set fill value to no_data_value
    nc_def_var_fill( ncid, varid_elevation, 0, &no_data_value );

    // end of header definition
    nc_enddef( ncid );

    // Write elevation data
    nc_put_var_double( ncid, varid_elevation, height_data_out.data() );

    // Write x data
    // Since we transpose the height data, we have to write our x_data to the netcdf "y" attribute
    // QGIS also expects the coordinates of pixel centers, so we shift by 0.5 * cell_size
    VectorX x_data_out = x_data + 0.5 * cell_size;
    nc_put_var_double( ncid, varid_y, x_data_out.data() );

    // Write y data
    // Since we transpose the height data and flip the y values, we have to write our y_data to the netcdf "x" attribute
    // and flip it
    // QGOS also expects the coordinates of pixel centers, so we shift by 0.5 * cell_size
    VectorX y_data_out = xt::flip( y_data ) + 0.5 * cell_size;
    nc_put_var_double( ncid, varid_x, y_data_out.data() );

    nc_close( ncid );
}

} // namespace Flowy
