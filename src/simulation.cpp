#include "simulation.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include "math.hpp"
#include "probability_dist.hpp"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#include "xtensor-blas/xlinalg.hpp"
#pragma GCC diagnostic pop
#include <algorithm>
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

void Simulation::run()
{
    for( int idx_flow = 0; idx_flow < input.n_flows; idx_flow++ )
    {
        // Determine n_lobes
        // @TODO: Calculate this later, using a beta distribution etc.
        int n_lobes = 0.5 * ( input.min_n_lobes + input.max_n_lobes );
        n_lobes_total += n_lobes;

        // Calculated for each flow with n_lobes number of lobes
        double delta_lobe_thickness
            = 2.0 * ( lobe_dimensions.avg_lobe_thickness - lobe_dimensions.thickness_min ) / ( n_lobes - 1.0 );

        // Build initial lobes which do not propagate descendents
        for( int idx_lobe = 0; idx_lobe < input.n_init; idx_lobe++ )
        {

            compute_initial_lobe_position( idx_flow, lobes[idx_lobe] );

            // Compute the thickness of the lobe
            lobes[idx_lobe].thickness = lobe_dimensions.thickness_min + idx_lobe * delta_lobe_thickness;

            auto [height_lobe_center, slope] = topography.height_and_slope( lobes[idx_lobe].center );

            // Perturb the angle (and set it)
            perturb_lobe_angle( lobes[idx_lobe], slope );
            // compute lobe axes
            compute_lobe_axes( lobes[idx_lobe], slope );

            // Add rasterized lobe
            topography.add_lobe( lobes[idx_lobe] );
        }
    }
}

} // namespace Flowtastic