// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/config.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/math.hpp"
#include "flowy/include/simulation.hpp"
#include <fmt/format.h>
#include <fmt/ostream.h>
#include <fmt/std.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_range_equals.hpp>
#include <cstddef>
#include <filesystem>
#include <vector>

TEST_CASE( "perturb_angle", "[perturb_angle]" )
{
    namespace fs = std::filesystem;
    using namespace Flowy;

    auto proj_root_path = fs::current_path();
    // Stupidly, we need an asc file to construct Simulation at the moment
    auto asc_file_path                  = proj_root_path / fs::path( "test/res/asc/file.asc" );
    auto input                          = Config::InputParams();
    input.source                        = asc_file_path;
    input.max_slope_prob                = 0.25;
    input.total_volume                  = 1.0;
    input.prescribed_avg_lobe_thickness = 1.0;

    INFO( fmt::format( "{}", fmt::streamed( input.source ) ) );
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
        my_lobe.set_azimuthal_angle( std::atan2( slope[1], slope[0] ) ); // Sets the angle prior to perturbation
        const double slope_norm = xt::linalg::norm( slope, 2 );          // Similar to np.linalg.norm
        simulation.perturb_lobe_angle( my_lobe, slope_norm );
        angle_samples[i] = my_lobe.get_azimuthal_angle();

        REQUIRE( ( ( my_lobe.get_azimuthal_angle() > -Math::pi ) && ( my_lobe.get_azimuthal_angle() < Math::pi ) ) );
    }

    double mean  = xt::mean( angle_samples )();
    double sigma = xt::stddev( angle_samples )();
    INFO( fmt::format( "\nmean = {},  mean_exp = {}\n", mean, mean_expected ) );
    INFO( fmt::format( "\nsigma = {},  sigma_gaussian = {}\n", sigma, sigma_gaussian ) );

    REQUIRE_THAT( mean, Catch::Matchers::WithinAbs( 0, 0.01 ) );
    REQUIRE( sigma < sigma_gaussian );
}

TEST_CASE( "budding_point", "[budding_point]" )
{
    using namespace Flowy;
    namespace fs = std::filesystem;

    Lobe lobe_parent;
    lobe_parent.center    = { -0.5, 0.5 };
    lobe_parent.semi_axes = { std::sqrt( 2 ) / 2, std::sqrt( 2 ) / 4 };
    lobe_parent.set_azimuthal_angle( -Math::pi / 4 );

    Lobe lobe_cur;
    lobe_cur.semi_axes = lobe_parent.semi_axes;

    auto proj_root_path = fs::current_path();
    auto asc_file_path  = proj_root_path / fs::path( "test/res/asc/file.asc" );
    Config::InputParams input_params;
    input_params.dist_fact                     = 1.0;
    input_params.source                        = asc_file_path;
    input_params.total_volume                  = 1;
    input_params.prescribed_avg_lobe_thickness = 1;

    Vector2 final_budding_point = lobe_parent.point_at_angle( 0 );

    auto simulation = Simulation( input_params, std::nullopt );
    simulation.compute_descendent_lobe_position( lobe_cur, lobe_parent, final_budding_point );

    INFO( fmt::format( "budding_point = {}\n", fmt::streamed( final_budding_point ) ) );
    INFO( fmt::format( "lobe_cur.center = {}\n", fmt::streamed( lobe_cur.center ) ) );

    Vector2 lobe_center_expected = { 0.5, -0.5 };

    REQUIRE( xt::isclose( lobe_cur.center, lobe_center_expected )() );
}

TEST_CASE( "hazard map", "[simulation]" )
{
    using namespace Flowy;
    namespace fs = std::filesystem;

    auto proj_root_path = fs::current_path();
    auto asc_file_path  = proj_root_path / fs::path( "test/res/asc/file.asc" );
    Config::InputParams input_params;
    input_params.dist_fact                     = 1.0;
    input_params.source                        = asc_file_path;
    input_params.total_volume                  = 1;
    input_params.prescribed_avg_lobe_thickness = 1;
    input_params.n_init                        = 1;

    auto lobes = std::vector<Lobe>( 10 );

    // Every lobe descends from the previous one
    for( size_t i = 1; i < lobes.size(); i++ )
    {
        lobes[i].idx_parent = i - 1;
    }

    // Except lobe 5 and 6
    lobes[5].idx_parent = 3;
    lobes[6].idx_parent = 4;

    auto simulation = Simulation( input_params, std::nullopt );
    simulation.compute_cumulative_descendents( lobes );

    std::vector<int> n_descendents_expected = { 9, 8, 7, 6, 4, 0, 3, 2, 1, 0 };

    for( size_t i = 0; i < lobes.size(); i++ )
    {
        INFO( fmt::format(
            "lobe[{}].n_descendents = {}, idx_parent = {}\n", i, lobes[i].n_descendents, lobes[i].idx_parent ) );
        REQUIRE( lobes[i].n_descendents == n_descendents_expected[i] );
    }
}
