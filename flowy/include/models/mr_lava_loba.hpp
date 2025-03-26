#pragma once

// GPL v3 License
// Copyright 2023--present Flowy developers, Mattia de' Michieli Vitturi, Simone Tarquini
#include "flowy/include/config.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/simulation.hpp"
#include "flowy/include/vent_flags.hpp"
#include "pdf_cpplib/include/probability_dist.hpp"
#include <algorithm>
#include <memory>
#include <vector>

namespace Flowy
{

class CommonLobeDimensions
{
public:
    CommonLobeDimensions() = default;
    explicit CommonLobeDimensions( const Config::InputParams & input )
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
        thickness_min     = 2.0 * input.thickness_ratio / ( input.thickness_ratio + 1.0 ) * avg_lobe_thickness;

        FLOWY_CHECK( lobe_area );
        FLOWY_CHECK( avg_lobe_thickness );
        FLOWY_CHECK( thickness_min );
        FLOWY_CHECK( exp_lobe_exponent );
        FLOWY_CHECK( max_semiaxis );
    }

    double avg_lobe_thickness = 0;
    double lobe_area          = 0;
    double max_semiaxis       = 0;
    double thickness_min      = 0;
    double exp_lobe_exponent  = 1;
};

class MrLavaLoba
{
public:
    CommonLobeDimensions lobe_dimensions;
    Config::InputParams input;

    MrLavaLoba( const Config::InputParams & input, std::mt19937 & gen )
            : lobe_dimensions( input ), input( input ), gen( gen )
    {
    }

    explicit MrLavaLoba( Simulation & simulation )
            : lobe_dimensions( simulation.input ), input( simulation.input ), gen( simulation.gen )
    {
    }

    // Calculate n_lobes
    int compute_n_lobes( int idx_flow ) const
    {
        int n_lobes{};
        // Number of lobes in the flow is a random number between the min and max values
        if( input.a_beta == 0 && input.b_beta == 0 )
        {
            std::uniform_int_distribution<> dist_num_lobes( input.min_n_lobes, input.max_n_lobes );
            n_lobes = dist_num_lobes( gen );
        }
        // Deterministic number of lobes, such that a beta probability density distribution is used (not a beta
        // distribution). However this means that n_lobes could potentially be greater than min_n_lobes
        else
        {
            double x_beta        = ( 1.0 * idx_flow ) / ( input.n_flows - 1.0 );
            double random_number = Math::beta_pdf( x_beta, input.a_beta, input.b_beta );
            n_lobes              = int(
                std::round( input.min_n_lobes + 0.5 * ( input.max_n_lobes - input.min_n_lobes ) * random_number ) );
        }
        return n_lobes;
    }

    // Calculate the current thickness of a lobe
    double compute_current_lobe_thickness( int idx_lobe, int n_lobes ) const
    {
        double delta_lobe_thickness
            = 2.0 * ( lobe_dimensions.avg_lobe_thickness - lobe_dimensions.thickness_min ) / ( n_lobes - 1.0 );

        // Calculated for each flow with n_lobes number of lobes
        FLOWY_CHECK( delta_lobe_thickness );

        return ( 1.0 - input.thickening_parameter )
               * ( lobe_dimensions.thickness_min + idx_lobe * delta_lobe_thickness );
    }

    // calculates the initial lobe position
    void compute_initial_lobe_position( int idx_flow, Lobe & lobe ) const
    {
        std::unique_ptr<VentFlag> f{};

        // Initial lobes are on the vent and flow starts from the first vent, second vent and so on
        if( input.vent_flag == 0 )
        {
            f = std::make_unique<VentFlag0>( idx_flow, input.n_flows, input.vent_coordinates, gen );
        }
        else if( input.vent_flag == 1 )
        {
            f = std::make_unique<VentFlag1>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 2 )
        {
            f = std::make_unique<VentFlag2>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 3 )
        {
            f = std::make_unique<VentFlag3>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 4 )
        {
            f = std::make_unique<VentFlag4>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 5 )
        {
            f = std::make_unique<VentFlag5>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 6 )
        {
            f = std::make_unique<VentFlag6>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 7 )
        {
            f = std::make_unique<VentFlag7>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }
        else if( input.vent_flag == 8 )
        {
            f = std::make_unique<VentFlag8>(
                input.vent_coordinates, input.fissure_probabilities, input.fissure_end_coordinates, gen );
        }

        lobe.center = f->get_position();
    }

    // perturbes the initial azimuthal angle of the lobe, which is
    // computed from the terrain slope
    void compute_lobe_axes( Lobe & lobe, double slope ) const
    {
        // Factor for the lobe eccentricity
        double aspect_ratio = std::min( input.max_aspect_ratio, 1.0 + input.aspect_ratio_coeff * slope );

        // Compute the semi-axes of the lobe
        double semi_major_axis = std::sqrt( lobe_dimensions.lobe_area / Math::pi ) * std::sqrt( aspect_ratio );
        double semi_minor_axis = std::sqrt( lobe_dimensions.lobe_area / Math::pi ) / std::sqrt( aspect_ratio );
        // Set the semi-axes
        lobe.semi_axes = { semi_major_axis, semi_minor_axis };
    }

    void compute_descendent_lobe_position( Lobe & lobe, const Lobe & parent, const Vector2 & final_budding_point ) const
    {
        Vector2 direction_to_new_lobe
            = ( final_budding_point - parent.center ) / xt::linalg::norm( final_budding_point - parent.center );
        Vector2 new_lobe_center = final_budding_point + input.dist_fact * direction_to_new_lobe * lobe.semi_axes[0];
        lobe.center             = new_lobe_center;
    }

    void perturb_lobe_angle( Lobe & lobe, double slope ) const
    {
        const double slope_deg = 180.0 / Math::pi * std::atan( slope );

        if( input.max_slope_prob < 1 )
        {
            if( slope_deg > 0.0 && input.max_slope_prob > 0 )
            {
                // Since we use radians instead of degrees, max_slope_prob has to be rescaled accordingly
                const double sigma = ( 1.0 - input.max_slope_prob ) / input.max_slope_prob * ( 90.0 - slope_deg )
                                     / slope_deg * Math::pi / 180.0;

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

    // Select which lobe amongst the existing lobes will be the parent for the new descendent lobe
    int select_parent_lobe( int idx_descendant, std::vector<Lobe> & lobes ) const
    {
        Lobe & lobe_descendent = lobes[idx_descendant];

        int idx_parent{};

        // Generate from the last lobe
        if( input.lobe_exponent <= 0 )
        {
            idx_parent = idx_descendant - 1;
        }
        else if( input.lobe_exponent >= 1 ) // Draw from a uniform random distribution if exponent is 1
        {
            std::uniform_int_distribution<int> dist_int( 0, idx_descendant - 1 );
            idx_parent = dist_int( gen );
        }
        else
        {
            std::uniform_real_distribution<double> dist( 0, 1 );
            const double idx0 = dist( gen );
            const auto idx1   = std::pow( idx0, input.lobe_exponent );
            idx_parent        = idx_descendant * idx1;
        }

        // Update the lobe information
        lobe_descendent.idx_parent   = idx_parent;
        lobe_descendent.dist_n_lobes = lobes[idx_parent].dist_n_lobes + 1;
        lobe_descendent.parent_weight *= lobe_dimensions.exp_lobe_exponent;

        return idx_parent;
    }

    void add_inertial_contribution( Lobe & lobe, const Lobe & parent, double slope ) const
    {
        double cos_angle_parent = parent.get_cos_azimuthal_angle();
        double sin_angle_parent = parent.get_sin_azimuthal_angle();
        double cos_angle_lobe   = lobe.get_cos_azimuthal_angle();
        double sin_angle_lobe   = lobe.get_sin_azimuthal_angle();

        double alpha_inertial = 0.0;

        const double eta = input.inertial_exponent;
        if( eta > 0 )
        {
            alpha_inertial = std::pow( ( 1.0 - std::pow( 2.0 * std::atan( slope ) / Math::pi, eta ) ), ( 1.0 / eta ) );
        }

        const double x_avg = ( 1.0 - alpha_inertial ) * cos_angle_lobe + alpha_inertial * cos_angle_parent;
        const double y_avg = ( 1.0 - alpha_inertial ) * sin_angle_lobe + alpha_inertial * sin_angle_parent;

        lobe.set_azimuthal_angle( std::atan2( y_avg, x_avg ) );
    }

private:
    std::mt19937 & gen;
};

} // namespace Flowy
