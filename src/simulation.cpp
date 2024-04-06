#include "simulation.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include "math.hpp"
#include "probability_dist.hpp"
#include "reservoir_sampling.hpp"
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace Flowy
{

CommonLobeDimensions::CommonLobeDimensions( const Config::InputParams & input, const AscFile & asc_file )
{
    if( !input.total_volume.has_value() )
        throw std::runtime_error( "Total volume flag not set" );

    if( input.fixed_dimension_flag == 1 )
    {
        if( !input.prescribed_lobe_area.has_value() )
            throw std::runtime_error( "prescribed_lobe_area is not set" );

        lobe_area          = input.prescribed_lobe_area.value();
        avg_lobe_thickness = input.total_volume.value()
                             / ( input.n_flows * input.prescribed_lobe_area.value() * 0.5
                                 * ( input.min_n_lobes + input.max_n_lobes ) );
    }
    else
    {
        if( !input.prescribed_avg_lobe_thickness.has_value() )
            throw std::runtime_error( "prescribed_avg_lobe_thickness is not set" );

        avg_lobe_thickness = input.prescribed_avg_lobe_thickness.value();
        lobe_area          = input.total_volume.value()
                    / ( input.n_flows * input.prescribed_avg_lobe_thickness.value() * 0.5
                        * ( input.min_n_lobes + input.max_n_lobes ) );
    }

    exp_lobe_exponent = std::exp( input.lobe_exponent );
    max_semiaxis      = std::sqrt( lobe_area * input.max_aspect_ratio / Math::pi );
    max_cells         = std::ceil( 2.0 * max_semiaxis / asc_file.cell_size ) + 2;
    thickness_min     = 2.0 * input.thickness_ratio / ( input.thickness_ratio + 1.0 ) * avg_lobe_thickness;
}

void Simulation::compute_initial_lobe_position( int idx_flow, Lobe & lobe )
{
    // For now, we've just implemented vent_flag = 0
    int idx_vent = std::floor( idx_flow * input.n_vents() / input.n_flows );
    lobe.center  = input.vent_coordinates[idx_vent];
}

void Simulation::write_lobe_data_to_file( const std::vector<Lobe> & lobes, const std::filesystem::path & path )
{
    std::fstream file;
    file.open( path, std::fstream::in | std::fstream::out | std::fstream::trunc );

    if( !file.is_open() )
    {
        throw std::runtime_error( fmt::format( "Unable to create output asc file: '{}'", path.string() ) );
    }

    file << fmt::format( "azimuthal_angle,centerx,centery,major_axis,minor_axis,dist_n_lobes,parent_weight,"
                         "n_descendents,idx_parent,alpha_intertial,thickness,height_center,slopex,slopey\n" );

    for( const auto & lobe : lobes )
    {

        file << fmt::format( "{},", lobe.get_azimuthal_angle() );
        file << fmt::format( "{},", lobe.center( 0 ) );
        file << fmt::format( "{},", lobe.center( 1 ) );
        file << fmt::format( "{},", lobe.semi_axes( 0 ) );
        file << fmt::format( "{},", lobe.semi_axes( 1 ) );
        file << fmt::format( "{},", lobe.dist_n_lobes );
        file << fmt::format( "{},", lobe.parent_weight );
        file << fmt::format( "{},", lobe.n_descendents );
        file << fmt::format( "{},", lobe.idx_parent.value_or( -1 ) );
        file << fmt::format( "{},", lobe.alpha_inertial );
        file << fmt::format( "{},", lobe.thickness );

        auto const [height, slope] = topography.height_and_slope( lobe.center );
        file << fmt::format( "{},", height );
        file << fmt::format( "{},", slope[0] );
        file << fmt::format( "{}", slope[1] );
        file << "\n";
    }
    file.close();
}

void Simulation::perturb_lobe_angle( Lobe & lobe, const Vector2 & slope )
{
    lobe.set_azimuthal_angle( std::atan2( slope[1], slope[0] ) ); // Sets the angle prior to perturbation
    const double slope_norm = xt::linalg::norm( slope, 2 );       // Similar to np.linalg.norm
    const double slope_deg  = std::atan( slope_norm );

    if( input.max_slope_prob < 1 )
    {
        if( slope_deg > 0.0 && input.max_slope_prob > 0 )
        {
            const double sigma
                = ( 1.0 - input.max_slope_prob ) / input.max_slope_prob * ( Math::pi / 2.0 - slope_deg ) / slope_deg;

            ProbabilityDist::truncated_normal_distribution<double> dist_truncated( 0, sigma, -Math::pi, Math::pi );
            const double angle_perturbation = dist_truncated( gen );
            lobe.set_azimuthal_angle( lobe.get_azimuthal_angle() + angle_perturbation );
        }
        else
        {
            std::uniform_real_distribution<double> dist_uniform( -Math::pi / 2, Math::pi / 2 );
            const double angle_perturbation = dist_uniform( gen );
            lobe.set_azimuthal_angle( lobe.get_azimuthal_angle() + angle_perturbation );
        }
    }
}

void Simulation::compute_lobe_axes( Lobe & lobe, const Vector2 & slope ) const
{
    const double slope_norm = xt::linalg::norm( slope, 2 );

    // Factor for the lobe eccentricity
    double aspect_ratio = std::min( input.max_aspect_ratio, 1.0 + input.aspect_ratio_coeff * slope_norm );

    // Compute the semi-axes of the lobe
    double semi_major_axis = std::sqrt( lobe_dimensions.lobe_area / Math::pi ) * std::sqrt( aspect_ratio );
    double semi_minor_axis = std::sqrt( lobe_dimensions.lobe_area / Math::pi ) / std::sqrt( aspect_ratio );
    // Set the semi-axes
    lobe.semi_axes = { semi_major_axis, semi_minor_axis };
}

// Select which lobe amongst the existing lobes will be the parent for the new descendent lobe
int Simulation::select_parent_lobe( int idx_descendant )
{

    // Their implementation
    // std::uniform_real_distribution<double> dist( 0, 1 );
    // double idx0 = dist( gen );
    // auto idx1 = std::pow( idx0, input.lobe_exponent );
    // auto idx2 = idx_descendant * idx1;
    // int idx_parent = std::floor( idx2 );
    // idx_parent = idx_descendant - 1;

    int idx_parent{};

    // Generate from the last lobe
    if( input.lobe_exponent == 0 )
    {
        idx_parent = idx_descendant - 1;
    }
    else if( input.lobe_exponent == 1 ) // Draw from a uniform random distribution if exponent is 1
    {
        std::uniform_int_distribution<int> dist_int( 0, idx_descendant - 1 );
        idx_parent = dist_int( gen );
    }
    else
    {
        // Otherwise, draw from an exponential distribution
        auto weight_probability_callback
            = [&]( int idx_current ) -> double { return lobes[idx_current].parent_weight; };

        idx_parent = ProbabilityDist::ReservoirSampling::reservoir_sampling_A_ExpJ(
            1, idx_descendant, weight_probability_callback, gen )[0];
    }

    // Update the lobe information
    lobes[idx_descendant].idx_parent   = idx_parent;
    lobes[idx_descendant].dist_n_lobes = lobes[idx_parent].dist_n_lobes + 1;
    lobes[idx_descendant].parent_weight *= lobe_dimensions.exp_lobe_exponent;

    // Loop over all ancestors and increase n_descendents
    std::optional<int> idx_parent_current = idx_parent;
    while( idx_parent_current.has_value() )
    {
        Lobe & lobe_current = lobes[idx_parent_current.value()];
        lobe_current.n_descendents++;
        idx_parent_current = lobe_current.idx_parent;
    }

    return idx_parent;
}

void Simulation::add_inertial_contribution( Lobe & lobe, const Lobe & parent, const Vector2 & slope ) const
{
    const double slope_norm = xt::linalg::norm( slope, 2 );
    double cos_angle_parent = parent.get_cos_azimuthal_angle();
    double sin_angle_parent = parent.get_sin_azimuthal_angle();
    double cos_angle_lobe   = lobe.get_cos_azimuthal_angle();
    double sin_angle_lobe   = lobe.get_sin_azimuthal_angle();

    double alpha_inertial = 0.0;

    const double eta = input.inertial_exponent;
    if( eta > 0 )
    {
        alpha_inertial = std::pow( ( 1.0 - std::pow( 2.0 * std::atan( slope_norm ) / Math::pi, eta ) ), ( 1.0 / eta ) );
    }

    const double x_avg = ( 1.0 - alpha_inertial ) * cos_angle_parent + alpha_inertial * cos_angle_lobe;
    const double y_avg = ( 1.0 - alpha_inertial ) * sin_angle_parent + alpha_inertial * sin_angle_lobe;

    lobe.set_azimuthal_angle( std::atan2( y_avg, x_avg ) );
}

void Simulation::compute_descendent_lobe_position( Lobe & lobe, const Lobe & parent, Vector2 final_budding_point )
{
    Vector2 direction_to_new_lobe
        = ( final_budding_point - parent.center ) / xt::linalg::norm( final_budding_point - parent.center );
    Vector2 new_lobe_center = final_budding_point + direction_to_new_lobe * lobe.semi_axes[0];
    lobe.center             = new_lobe_center;
}

bool Simulation::stop_condition( const Vector2 & point, double radius )
{
    return topography.is_point_near_boundary( point, radius )
           || topography.get_height( point ) < asc_file.no_data_value + 1;
}

void Simulation::run()
{
    int n_lobes_total     = 0; // This is the total number of lobes, accumulated over all flows
    int n_lobes_processed = 0;

    // Save initial topography to asc file
    auto asc_file = topography.to_asc_file();
    asc_file.save( input.output_folder / "initial.asc" );

    auto t_run_start = std::chrono::high_resolution_clock::now();
    for( int idx_flow = 0; idx_flow < input.n_flows; idx_flow++ )
    {
        // Determine n_lobes
        // @TODO: Calculate this later, using a beta distribution etc.
        int n_lobes = 0.5 * ( input.min_n_lobes + input.max_n_lobes );
        n_lobes_total += n_lobes;

        lobes = std::vector<Lobe>{};
        lobes.reserve( n_lobes );

        // Calculated for each flow with n_lobes number of lobes
        double delta_lobe_thickness
            = 2.0 * ( lobe_dimensions.avg_lobe_thickness - lobe_dimensions.thickness_min ) / ( n_lobes - 1.0 );

        // Build initial lobes which do not propagate descendents
        for( int idx_lobe = 0; idx_lobe < input.n_init; idx_lobe++ )
        {
            lobes.emplace_back();
            Lobe & lobe_cur = lobes[idx_lobe];

            compute_initial_lobe_position( idx_flow, lobe_cur );

            // Compute the thickness of the lobe
            lobe_cur.thickness = lobe_dimensions.thickness_min + idx_lobe * delta_lobe_thickness;

            auto [height_lobe_center, slope] = topography.height_and_slope( lobe_cur.center );

            // Perturb the angle (and set it)
            perturb_lobe_angle( lobe_cur, slope );

            // compute lobe axes
            compute_lobe_axes( lobe_cur, slope );

            // Add rasterized lobe
            topography.add_lobe( lobe_cur );

            n_lobes_processed++;
        }

        // Loop over the rest of the lobes (skipping the initial ones).
        // Each lobe is a descendant of a parent lobe
        for( int idx_lobe = input.n_init; idx_lobe < n_lobes; idx_lobe++ )
        {
            lobes.emplace_back();
            Lobe & lobe_cur = lobes[idx_lobe];

            // Select which of the previously created lobes is the parent lobe
            // from which the new descendent lobe will bud
            auto idx_parent = select_parent_lobe( idx_lobe );

            Lobe & lobe_parent = lobes[idx_parent];

            // stopping condition (parent lobe close the domain boundary or at a not defined z value)
            if( stop_condition( lobe_parent.center, lobe_parent.semi_axes[0] ) )
            {
                break;
            }

            // Find the preliminary budding point on the perimeter of the parent lobe (npoints is the number of raster
            // points on the ellipse)
            Flowy::Vector2 budding_point = topography.find_preliminary_budding_point( lobe_parent, input.npoints );

            auto [height_lobe_center, slope_parent] = topography.height_and_slope( lobe_parent.center );

            // Perturb the angle and set it (not on the parent anymore)
            perturb_lobe_angle( lobe_cur, slope_parent );

            // Add the inertial contribution
            add_inertial_contribution( lobe_cur, lobe_parent, slope_parent );

            // Compute the final budding point
            // It is defined by the point on the perimeter of the parent lobe closest to the center of the new lobe
            Vector2 final_budding_point
                = 2.0 * lobe_parent.center - lobe_parent.point_at_angle( lobe_cur.get_azimuthal_angle() );
            // Vector2 final_budding_point = lobe_parent.point_at_angle( lobe_cur.get_azimuthal_angle() );

            if( stop_condition( final_budding_point ) )
            {
                break;
            }
            // Get the slope at the final budding point
            auto [height_budding_point, slope_budding_point] = topography.height_and_slope( final_budding_point );

            // compute the new lobe axes
            compute_lobe_axes( lobe_cur, slope_budding_point );

            // Get new lobe center
            compute_descendent_lobe_position( lobe_cur, lobe_parent, final_budding_point );

            if( stop_condition( lobe_cur.center ) )
            {
                break;
            }

            // Compute the thickness of the lobe
            lobe_cur.thickness = lobe_dimensions.thickness_min + idx_lobe * delta_lobe_thickness;

            // Add rasterized lobe
            topography.add_lobe( lobe_cur );
            n_lobes_processed++;
        }

        // write_lobe_data_to_file( lobes, fmt::format( "flow_{}_lobes_{}.txt", idx_flow, n_lobes_processed ) );

        auto t_cur          = std::chrono::high_resolution_clock::now();
        auto remaining_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            ( input.n_flows - idx_flow - 1 ) * ( t_cur - t_run_start ) / ( idx_flow + 1 ) );
        fmt::print( "     remaining_time = {:%Hh %Mm %Ss}\n", remaining_time );
    }

    auto t_cur      = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>( ( t_cur - t_run_start ) );
    fmt::print( "total_time = {:%Hh %Mm %Ss}\n", total_time );

    fmt::print( "Total number of processed lobes = {}\n", n_lobes_processed );

    fmt::print( "n_lobes/ms = {}\n", ( n_lobes_processed / total_time.count() ) );

    // Save final topography to asc file
    asc_file = topography.to_asc_file();
    asc_file.save( input.output_folder / "output.asc" );
}

} // namespace Flowy