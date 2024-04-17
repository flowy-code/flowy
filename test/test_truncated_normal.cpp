#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include "math.hpp"
#include "pdf_cpplib/include/probability_dist.hpp"
#include <fmt/ostream.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
namespace fs = std::filesystem;

template<typename T, std::size_t N>
std::ostream & operator<<( std::ostream & os, std::array<T, N> const & v1 )
{
    std::for_each( begin( v1 ), end( v1 ), [&os]( T val ) { os << val << " "; } );
    return os;
}

// Samples the distribution n_samples times and writes results to file
template<typename distT>
void write_results_to_file( int N_Samples, distT dist, const std::string & filename )
{
    auto proj_root_path = fs::current_path();
    auto file           = proj_root_path / fs::path( "test/output_probability_distributions/" + filename );
    fmt::print( "file = {}\n", fmt::streamed( file ) );
    fs::create_directories( file.parent_path() );

    auto gen = std::mt19937( 0 );

    std::ofstream filestream( file );
    filestream << std::setprecision( 16 );

    for( int i = 0; i < N_Samples; i++ )
    {
        filestream << dist( gen ) << "\n";
    }
    filestream.close();
}

TEST_CASE( "Test the probability distributions", "[prob]" )
{
    write_results_to_file(
        10000, ProbabilityDist::truncated_normal_distribution<double>( 0.0, 2, -Flowy::Math::pi, Flowy::Math::pi ),
        "truncated_normal.txt" );
}