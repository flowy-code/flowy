// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/asc_file.hpp"
#include "flowy/include/dump_csv.hpp"
#include "flowy/include/topography_file.hpp"
#include <fast_float/fast_float.h>
#include <fmt/compile.h>
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

    // Fast bulk parse with a locale-independent from_chars-style API.
    std::string buf( ( std::istreambuf_iterator<char>( file ) ), std::istreambuf_iterator<char>() );
    std::vector<double> vals;
    vals.reserve( nrows_header * ncols_header );
    const char * p   = buf.data();
    const char * end = p + buf.size();
    while( p < end )
    {
        while( p < end && ( *p == ' ' || *p == '\n' || *p == '\r' || *p == '\t' ) )
            ++p;
        if( p >= end )
            break;
        double d{};
        auto [next, ec] = fast_float::from_chars( p, end, d );
        if( ec != std::errc() )
        {
            ++p;
            continue;
        }
        vals.push_back( d );
        p = next;
    }
    std::array<std::size_t, 2> shp = { nrows_header, ncols_header };
    data                           = xt::adapt( vals, shp );

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

    // Bulk format into one buffer, then write it once.
    const size_t ncols = data.shape()[0];
    const size_t nrows = data.shape()[1];

    std::string buf;
    buf.reserve( ncols * nrows * 15 + 256 );
    auto out = std::back_inserter( buf );
    out      = fmt::format_to( out, FMT_COMPILE( "ncols {}\n" ), ncols );
    out      = fmt::format_to( out, FMT_COMPILE( "nrows {}\n" ), nrows );
    out      = fmt::format_to( out, FMT_COMPILE( "xllcorner {}\n" ), lower_left_corner()[0] );
    out      = fmt::format_to( out, FMT_COMPILE( "yllcorner {}\n" ), lower_left_corner()[1] );
    out      = fmt::format_to( out, FMT_COMPILE( "cellsize {}\n" ), cell_size() );
    out      = fmt::format_to( out, FMT_COMPILE( "NODATA_value {}\n" ), no_data_value );

    for( size_t r = 0; r < nrows; r++ )
    {
        const size_t src_row = nrows - 1 - r; // undo the y-flip
        for( size_t c = 0; c < ncols; c++ )
        {
            if( c > 0 )
                buf.push_back( ' ' );
            const double value = data( c, src_row );
            if( value == 0.0 )
                buf.push_back( '0' );
            else
                out = fmt::format_to( out, FMT_COMPILE( "{}" ), value );
        }
        buf.push_back( '\n' );
    }

    std::ofstream out_file( path, std::ios::binary );
    if( !out_file.is_open() )
        throw std::runtime_error( fmt::format( "Unable to create output asc file: '{}'", path.string() ) );
    out_file.write( buf.data(), buf.size() );
}

} // namespace Flowy
