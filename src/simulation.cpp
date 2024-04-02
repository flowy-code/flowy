#include "simulation.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include "math.hpp"
#include "probability_dist.hpp"
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <stdexcept>

namespace Flowtastic
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

    max_semiaxis  = std::sqrt( lobe_area * input.max_aspect_ratio / Math::pi );
    max_cells     = std::ceil( 2.0 * max_semiaxis / asc_file.cell_size ) + 2;
    thickness_min = 2.0 * input.thickness_ratio / ( input.thickness_ratio + 1.0 ) * avg_lobe_thickness;
}

void Simulation::compute_initial_lobe_position( int idx_flow, Lobe & lobe )
{
    // For now, we've just implemented vent_flag = 0
    int idx_vent = std::floor( idx_flow * input.n_vents() / input.n_flows );
    lobe.center  = input.vent_coordinates[idx_vent];
}

void Simulation::perturb_lobe_angle( Lobe & lobe, const Vector2 & slope )
{
    lobe.azimuthal_angle    = std::atan2( slope[1], slope[0] ); // Sets the angle prior to perturbation
    const double slope_norm = xt::linalg::norm( slope, 2 );     // Similar to np.linalg.norm
    const double slope_deg  = std::atan( slope_norm );

    if( input.max_slope_prob < 1 )
    {
        if( slope_deg > 0.0 && input.max_slope_prob > 0 )
        {
            const double sigma
                = ( 1.0 - input.max_slope_prob ) / input.max_slope_prob * ( Math::pi / 2.0 - slope_deg ) / slope_deg;

            ProbabilityDist::truncated_normal_distribution<double> dist_truncated( 0, sigma, -Math::pi, Math::pi );
            const double angle_perturbation = dist_truncated( gen );
            lobe.azimuthal_angle += angle_perturbation;
        }
        else
        {
            std::uniform_real_distribution<double> dist_uniform( -Math::pi, Math::pi );
            const double angle_perturbation = dist_uniform( gen );
            lobe.azimuthal_angle += angle_perturbation;
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
    // TODO: implement input.lobe_exponent, input.force_max_length and input.start_from_dist_flag
    std::uniform_int_distribution<int> dist_int( 0, idx_descendant - 1 );
    int idx_parent = dist_int( gen );

    lobes[idx_descendant].idx_parent   = idx_parent;
    lobes[idx_descendant].dist_n_lobes = lobes[idx_parent].dist_n_lobes + 1;

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
    double cos_angle_parent = std::cos( parent.azimuthal_angle );
    double sin_angle_parent = std::sin( parent.azimuthal_angle );
    double cos_angle_lobe   = std::cos( lobe.azimuthal_angle );
    double sin_angle_lobe   = std::sin( lobe.azimuthal_angle );

    double alpha_inertial = 0.0;

    const double eta = input.inertial_exponent;
    if( eta > 0 )
    {
        alpha_inertial = std::pow( ( 1.0 - std::pow( 2.0 * std::atan( slope_norm ) / Math::pi, eta ) ), ( 1.0 / eta ) );
    }

    const double x_avg = ( 1.0 - alpha_inertial ) * cos_angle_parent + alpha_inertial * cos_angle_lobe;
    const double y_avg = ( 1.0 - alpha_inertial ) * sin_angle_parent + alpha_inertial * sin_angle_lobe;

    lobe.azimuthal_angle = std::atan2( y_avg, x_avg );
}

void Simulation::compute_descendent_lobe_position( Lobe & lobe, const Lobe & parent, Vector2 final_budding_point )
{
    Vector2 direction_to_new_lobe
        = ( final_budding_point - parent.center ) / xt::linalg::norm( final_budding_point - parent.center );
    Vector2 new_lobe_center = final_budding_point + direction_to_new_lobe * lobe.semi_axes[0];
    lobe.center             = new_lobe_center;
}

void Simulation::run()
{
    for( int idx_flow = 0; idx_flow < input.n_flows; idx_flow++ )
    {
        auto t_flow_start = std::chrono::high_resolution_clock::now();

        // Determine n_lobes
        // @TODO: Calculate this later, using a beta distribution etc.
        int n_lobes = 0.5 * ( input.min_n_lobes + input.max_n_lobes );
        n_lobes_total += n_lobes;

        lobes = std::vector<Lobe>( n_lobes );

        // Calculated for each flow with n_lobes number of lobes
        double delta_lobe_thickness
            = 2.0 * ( lobe_dimensions.avg_lobe_thickness - lobe_dimensions.thickness_min ) / ( n_lobes - 1.0 );

        // Build initial lobes which do not propagate descendents
        for( int idx_lobe = 0; idx_lobe < input.n_init; idx_lobe++ )
        {
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
        }

        // Loop over the rest of the lobes (skipping the initial ones).
        // Each lobe is a descendant of a parent lobe
        for( int idx_lobe = input.n_init; idx_lobe < n_lobes; idx_lobe++ )
        {
            Lobe & lobe_cur = lobes[idx_lobe];

            // Select which of the previously created lobes is the parent lobe
            // from which the new descendent lobe will bud
            auto idx_parent = select_parent_lobe( idx_lobe );

            Lobe & lobe_parent = lobes[idx_parent];

            // Find the preliminary budding point on the perimeter of the parent lobe (npoints is the number of raster
            // points on the ellipse)
            Flowtastic::Vector2 budding_point = topography.find_preliminary_budding_point( lobe_parent, input.npoints );

            auto [height_lobe_center, slope_parent] = topography.height_and_slope( lobe_parent.center );

            // Perturb the angle and set it (not on the parent anymore)
            perturb_lobe_angle( lobe_cur, slope_parent );

            // Add the inertial contribution
            add_inertial_contribution( lobe_cur, lobe_parent, slope_parent );

            // Compute the final budding point
            // It is defined by the point on the perimeter of the parent lobe closest to the center of the new lobe
            Vector2 final_budding_point = lobe_parent.point_at_angle( lobe_cur.azimuthal_angle );

            // Get the slope at the final budding point
            auto [height_budding_point, slope_budding_point] = topography.height_and_slope( final_budding_point );

            // compute the new lobe axes
            compute_lobe_axes( lobe_cur, slope_budding_point );

            // Get new lobe center
            compute_descendent_lobe_position( lobe_cur, lobe_parent, final_budding_point );

            // Compute the thickness of the lobe
            lobe_cur.thickness = lobe_dimensions.thickness_min + idx_lobe * delta_lobe_thickness;

            // Add rasterized lobe
            topography.add_lobe( lobe_cur );
        }

        auto t_flow_end     = std::chrono::high_resolution_clock::now();
        auto remaining_time = std::chrono::duration_cast<std::chrono::seconds>(
            ( input.n_flows - idx_flow ) * ( t_flow_end - t_flow_start ) );

        fmt::print( "remaining_time = {:%Hh %Mm %Ss}\n", remaining_time );
    }
}

} // namespace Flowtastic