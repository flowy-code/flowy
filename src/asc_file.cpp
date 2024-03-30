#include "asc_file.hpp"
#include <fmt/format.h>
#include <algorithm>
#include <fstream>
#include <istream>
#include <stdexcept>
#include <xtensor/xcsv.hpp>

namespace Flowtastic
{
AscFile::AscFile( const std::filesystem::path & path )
{

    std::ifstream file( path.string() ); // Open the file
    if( !file.is_open() )
    {
        throw std::runtime_error( "Unable to open file" );
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

    auto ncols_header = std::stoi( get_number_string() );
    auto nrows_header = std::stoi( get_number_string() );

    auto lx           = std::stod( get_number_string() );
    auto ly           = std::stod( get_number_string() );
    lower_left_corner = { lx, ly };

    cell_size     = std::stod( get_number_string() );
    no_data_value = std::stod( get_number_string() );

    height_data = xt::load_csv<double>( file, ' ' );

    if( nrows_header != nrows() )
    {
        throw std::runtime_error(
            fmt::format( "nrows in header is {}, but there are {} rows of data", nrows_header, nrows() ) );
    }

    if( ncols_header != ncols() )
    {
        throw std::runtime_error(
            fmt::format( "ncols in header is {}, but there are {} cols of data", ncols_header, ncols() ) );
    }
}
} // namespace Flowtastic