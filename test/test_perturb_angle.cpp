#include "config.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include "math.hpp"
#include "simulation.hpp"
#include "xtensor-blas/xlinalg.hpp"
#include "xtensor/xbuilder.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xtensor_forward.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <filesystem>

TEST_CASE( "perturb_angle", "[perturb_angle]" )
{
    namespace fs = std::filesystem;
    using namespace Flowtastic;

    auto proj_root_path = fs::current_path();
    // Stupidly, we need an asc file to construct Simulation at the moment
    auto asc_file_path                  = proj_root_path / fs::path( "test/res/asc/file.asc" );
    auto input                          = Config::InputParams();
    input.source                        = asc_file_path;
    input.max_slope_prob                = 0.25;
    input.total_volume                  = 1.0;
    input.prescribed_avg_lobe_thickness = 1.0;

    fmt::print( "{}", fmt::streamed( input.source ) );
    auto simulation = Simulation( input, 0 );

    Lobe my_lobe;

    int n_samples                    = 100000;
    xt::xarray<double> angle_samples = xt::empty<double>( { n_samples } );

    Vector2 slope              = { 1.5, 0 };
    const double mean_expected = 0.0;

    const double slope_norm = xt::linalg::norm( slope, 2 ); // Similar to np.linalg.norm
    const double slope_deg  = std::atan( slope_norm );

    // NOTE: this is not the sigma we expect from the samples, but the sigma of the gaussian *before* truncation
    const double sigma_gaussian
        = ( 1.0 - input.max_slope_prob ) / input.max_slope_prob * ( Math::pi / 2.0 - slope_deg ) / slope_deg;

    for( int i = 0; i < n_samples; i++ )
    {
        simulation.perturb_lobe_angle( my_lobe, slope );
        angle_samples[i] = my_lobe.azimuthal_angle;

        REQUIRE( ( ( my_lobe.azimuthal_angle > -Math::pi ) && ( my_lobe.azimuthal_angle < Math::pi ) ) );
    }

    double mean  = xt::mean( angle_samples )();
    double sigma = xt::stddev( angle_samples )();
    fmt::print( "\nmean = {},  mean_exp = {}\n", mean, mean_expected );
    fmt::print( "\nsigma = {},  sigma_gaussian = {}\n", sigma, sigma_gaussian );

    REQUIRE_THAT( mean, Catch::Matchers::WithinAbs( 0, 0.01 ) );
    REQUIRE( sigma < sigma_gaussian );
}