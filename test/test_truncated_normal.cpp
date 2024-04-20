// GPL v3 License
// Copyright 2023--present Flowy developers
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>

#include "flowy/include/definitions.hpp"
#include "flowy/include/math.hpp"
#include "pdf_cpplib/include/probability_dist.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include <fmt/ostream.h>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
namespace fs = std::filesystem;

// Samples the distribution n_samples times and writes results to file
auto write_results_to_file( const Flowy::VectorX & samples, const std::string & filename )
{
    auto proj_root_path = fs::current_path();
    auto file           = proj_root_path / fs::path( "test/output_probability_distributions/" + filename );
    INFO( fmt::format( "file = {}\n", fmt::streamed( file ) ) );
    fs::create_directories( file.parent_path() );

    std::ofstream filestream( file );
    filestream << std::setprecision( 16 );

    for( size_t i = 0; i < samples.size(); i++ )
    {
        filestream << samples[i] << "\n";
    }
    filestream.close();
    return samples;
}

TEST_CASE( "Test the probability distributions", "[prob]" )
{
    int N_Samples = 100000;
    auto gen      = std::mt19937( 0 );
    auto dist = ProbabilityDist::truncated_normal_distribution<double>( 0.0, 2.0, -Flowy::Math::pi, Flowy::Math::pi );
    Flowy::VectorX samples = xt::empty<double>( { N_Samples } );

    for( size_t i = 0; i < samples.size(); i++ )
    {
        samples[i] = dist( gen );
    }

    // Use for debugging
    // write_results_to_file( samples, "truncated_normal.txt" );

    REQUIRE_THAT( xt::mean( samples )(), Catch::Matchers::WithinAbs( 0, 0.01 ) );
    REQUIRE_THAT( xt::amax( samples )(), Catch::Matchers::WithinAbs( Flowy::Math::pi, 0.01 ) );
    REQUIRE_THAT( xt::amin( samples )(), Catch::Matchers::WithinAbs( -Flowy::Math::pi, 0.01 ) );
}
