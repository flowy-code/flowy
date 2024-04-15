// GPL v3 License
// Copyright 2023--present Flowy developers
#include "include/asc_file.hpp"
#include "include/dump_csv.hpp"
#include <fmt/format.h>
#include <fstream>

namespace Flowy
{

AscFile::AscFile( const std::filesystem::path & path, std::optional<AscCrop> crop )
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

    auto lx           = std::stod( get_number_string() );
    auto ly           = std::stod( get_number_string() );
    lower_left_corner = { lx, ly };

    cell_size     = std::stod( get_number_string() );
    no_data_value = std::stod( get_number_string() );

    height_data = xt::load_csv<double>( file, ' ' );

    if( nrows_header != height_data.shape()[0] )
    {
        throw std::runtime_error( fmt::format(
            "nrows in header is {}, but there are {} rows of data", nrows_header, height_data.shape()[0] ) );
    }

    if( ncols_header != height_data.shape()[1] )
    {
        throw std::runtime_error( fmt::format(
            "ncols in header is {}, but there are {} cols of data", ncols_header, height_data.shape()[1] ) );
    }

    // Now we transform the height data to the shape that the rest of the code expects
    // That means first, we flip the direction of the y-axis and then we transpose
    height_data = xt::transpose( xt::flip( height_data, 0 ) );

    // If cropping is used, we slice the height_data array
    if( crop.has_value() )
    {
        int idx_x_min = std::clamp<int>( ( crop->x_min - lx ) / cell_size, 0, height_data.shape()[0] - 1 );
        int idx_x_max = std::clamp<int>( ( crop->x_max - lx ) / cell_size, 0, height_data.shape()[0] - 1 );
        int idx_y_min = std::clamp<int>( ( crop->y_min - ly ) / cell_size, 0, height_data.shape()[1] - 1 );
        int idx_y_max = std::clamp<int>( ( crop->y_max - ly ) / cell_size, 0, height_data.shape()[1] - 1 );

        height_data
            = xt::view( height_data, xt::range( idx_x_min, idx_x_max + 1 ), xt::range( idx_y_min, idx_y_max + 1 ) );

        lower_left_corner = { lx + idx_x_min * cell_size, ly + idx_y_min * cell_size };
    }

    this->x_data = xt::arange(
        lower_left_corner[0], lower_left_corner[0] + ( static_cast<double>( height_data.shape()[0] ) ) * cell_size,
        cell_size );

    this->y_data = xt::arange(
        lower_left_corner[1], lower_left_corner[1] + ( static_cast<double>( height_data.shape()[1] ) ) * cell_size,
        cell_size );
}

void AscFile::save( const std::filesystem::path & path )
{
    std::fstream file;
    file.open( path, std::fstream::in | std::fstream::out | std::fstream::trunc );

    if( !file.is_open() )
    {
        throw std::runtime_error( fmt::format( "Unable to create output asc file: '{}'", path.string() ) );
    }

    file << fmt::format( "ncols {}\n", height_data.shape()[0] );
    file << fmt::format( "nrows {}\n", height_data.shape()[1] );
    file << fmt::format( "xllcorner {}\n", lower_left_corner[0] );
    file << fmt::format( "yllcorner {}\n", lower_left_corner[1] );
    file << fmt::format( "cellsize {}\n", cell_size );
    file << fmt::format( "NODATA_value {}\n", no_data_value );

    // We have to undo the transformation we applied to the height data
    // Therefore, we transpose and *then* we flip the y-axis
    auto height_data_out = xt::flip( xt::transpose( height_data ), 0 );
    Utility::dump_csv( file, height_data_out, ' ' );

    file.close();
}

} // namespace Flowy
