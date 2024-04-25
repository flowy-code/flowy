// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/asc_file.hpp"
#include "flowy/include/dump_csv.hpp"
#include "flowy/include/topography_file.hpp"
#include <fmt/format.h>
#include <fstream>

namespace Flowy
{

AscFile::AscFile( const std::filesystem::path & path, const std::optional<TopographyCrop> & crop )
{
    std::ifstream file( path.string() ); // Open the file
    if( !file.is_open() )
    {
        throw std::runtime_error( fmt::format( "Unable to open asc file: '{}'", path.string() ) );
    }

    /* This is what the first six lines look like
    ncols 2
    nrows 2
    xllcorner 2.701332e+05
    yllcorner 2.123588e+06
    cellsize 20
    NODATA_value -9999
    */

    std::string line;
    auto get_number_string = [&]()
    {
        std::getline( file, line );
        auto pos_space = line.find( ' ', 0 );
        return line.substr( pos_space, std::string::npos );
    };

    size_t ncols_header = std::stoi( get_number_string() );
    size_t nrows_header = std::stoi( get_number_string() );

    // Parse lower left corner
    double lx = std::stod( get_number_string() );
    double ly = std::stod( get_number_string() );

    // Parse cell size
    double cell_size_file = std::stod( get_number_string() );

    no_data_value = std::stod( get_number_string() );

    data = xt::load_csv<double>( file, ' ' );

    if( nrows_header != data.shape()[0] )
    {
        throw std::runtime_error(
            fmt::format( "nrows in header is {}, but there are {} rows of data", nrows_header, data.shape()[0] ) );
    }

    if( ncols_header != data.shape()[1] )
    {
        throw std::runtime_error(
            fmt::format( "ncols in header is {}, but there are {} cols of data", ncols_header, data.shape()[1] ) );
    }

    // Now we transform the height data to the shape that the rest of the code expects
    // That means first, we flip the direction of the y-axis and then we transpose
    data = xt::transpose( xt::flip( data, 0 ) );

    this->x_data = xt::arange( lx, lx + ( static_cast<double>( data.shape()[0] ) ) * cell_size_file, cell_size_file );

    this->y_data = xt::arange( ly, ly + ( static_cast<double>( data.shape()[1] ) ) * cell_size_file, cell_size_file );

    if( crop.has_value() )
    {
        crop_topography( crop.value() );
    }
}

void AscFile::save( const std::filesystem::path & path_ )
{
    auto path = handle_suffix( path_ );

    std::fstream file;
    file.open( path, std::fstream::in | std::fstream::out | std::fstream::trunc );

    if( !file.is_open() )
    {
        throw std::runtime_error( fmt::format( "Unable to create output asc file: '{}'", path.string() ) );
    }

    file << fmt::format( "ncols {}\n", data.shape()[0] );
    file << fmt::format( "nrows {}\n", data.shape()[1] );
    file << fmt::format( "xllcorner {}\n", lower_left_corner()[0] );
    file << fmt::format( "yllcorner {}\n", lower_left_corner()[1] );
    file << fmt::format( "cellsize {}\n", cell_size() );
    file << fmt::format( "NODATA_value {}\n", no_data_value );

    // We have to undo the transformation we applied to the height data
    // Therefore, we transpose and *then* we flip the y-axis
    auto data_out = xt::flip( xt::transpose( data ), 0 );
    Utility::dump_csv( file, data_out, ' ' );

    file.close();
}

} // namespace Flowy
