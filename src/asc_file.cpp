#include <charconv>
#include <fcntl.h>
#include <unistd.h>
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

    // Fast bulk parse: read the rest of the file and from_chars it (replaces the
    // strtod/stream path that dominated the profile).
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
        auto [next, ec] = std::from_chars( p, end, d );
        if( ec != std::errc() )
        {
            ++p;
            continue;
        }
        vals.push_back( d );
        p = next;
    }
    std::array<std::size_t, 2> shp = { nrows_header, ncols_header };
    data = xt::adapt( vals, shp );

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

    // Bulk format into one buffer, then a single write() (replaces the per-value
    // stream formatting that dominated the write profile).
    const size_t ncols = data.shape()[0];
    const size_t nrows = data.shape()[1];

    std::string buf;
    buf.reserve( ncols * nrows * 15 + 256 );
    auto out = std::back_inserter( buf );
    out = fmt::format_to( out, "ncols {}\n", ncols );
    out = fmt::format_to( out, "nrows {}\n", nrows );
    out = fmt::format_to( out, "xllcorner {}\n", lower_left_corner()[0] );
    out = fmt::format_to( out, "yllcorner {}\n", lower_left_corner()[1] );
    out = fmt::format_to( out, "cellsize {}\n", cell_size() );
    out = fmt::format_to( out, "NODATA_value {}\n", no_data_value );

    for( size_t r = 0; r < nrows; r++ )
    {
        const size_t src_row = nrows - 1 - r; // undo the y-flip
        for( size_t c = 0; c < ncols; c++ )
        {
            if( c > 0 )
                buf.push_back( ' ' );
            out = fmt::format_to( out, "{}", data( c, src_row ) );
        }
        buf.push_back( '\n' );
    }

    int fd = ::open( path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644 );
    if( fd < 0 )
        throw std::runtime_error( fmt::format( "Unable to create output asc file: '{}'", path.string() ) );
    const char * wp = buf.data();
    size_t rem      = buf.size();
    while( rem > 0 )
    {
        ssize_t w = ::write( fd, wp, rem );
        if( w < 0 ) { ::close( fd ); throw std::runtime_error( "write error" ); }
        wp += w; rem -= static_cast<size_t>( w );
    }
    ::close( fd );
}

} // namespace Flowy
