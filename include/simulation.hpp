#pragma once
#include "asc_file.hpp"
#include "config.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include "topography.hpp"
#include <filesystem>
#include <random>
#include <vector>

namespace Flowy
{

class CommonLobeDimensions
{
public:
    CommonLobeDimensions() = default;
    CommonLobeDimensions( const Config::InputParams & input, const AscFile & asc_file );

    double avg_lobe_thickness = 0;
    double lobe_area          = 0;
    double max_semiaxis       = 0;
    int max_cells             = 0;
    double thickness_min      = 0;
    double exp_lobe_exponent  = 1;
};

class Simulation
{
public:
    Simulation( const Config::InputParams & input, std::optional<int> rng_seed ) : input( input )
    {

        this->rng_seed = rng_seed.value_or( std::random_device()() );
        gen            = std::mt19937( this->rng_seed );

        // Create output directory
        std::filesystem::create_directories( input.output_folder ); // Create the output directory

        // Crop if all of these have a value
        if( input.east_to_vent.has_value() && input.west_to_vent.has_value() && input.south_to_vent.has_value()
            && input.north_to_vent.has_value() )
        {
            auto crop = AscCrop{};

            auto min_x_it = std::min_element(
                input.vent_coordinates.begin(), input.vent_coordinates.end(),
                [&]( const Vector2 & p1, const Vector2 & p2 ) { return p1[0] < p2[0]; } );
            auto min_y_it = std::min_element(
                input.vent_coordinates.begin(), input.vent_coordinates.end(),
                [&]( const Vector2 & p1, const Vector2 & p2 ) { return p1[1] < p2[1]; } );
            auto max_x_it = std::max_element(
                input.vent_coordinates.begin(), input.vent_coordinates.end(),
                [&]( const Vector2 & p1, const Vector2 & p2 ) { return p1[0] < p2[0]; } );
            auto max_y_it = std::max_element(
                input.vent_coordinates.begin(), input.vent_coordinates.end(),
                [&]( const Vector2 & p1, const Vector2 & p2 ) { return p1[1] < p2[1]; } );

            crop.x_min = ( *min_x_it )[0] - input.west_to_vent.value();
            crop.x_max = ( *max_x_it )[0] + input.east_to_vent.value();
            crop.y_min = ( *min_y_it )[1] - input.south_to_vent.value();
            crop.y_max = ( *max_y_it )[1] + input.north_to_vent.value();

            asc_file = AscFile( input.source, crop );
        }
        else
        {
            asc_file = AscFile( input.source );
        }

        topography      = Topography( asc_file );
        lobe_dimensions = CommonLobeDimensions( input, asc_file );
    };

    Config::InputParams input;
    AscFile asc_file;
    Topography topography;
    CommonLobeDimensions lobe_dimensions;

    std::vector<Lobe> lobes; // Lobes per flows

    // calculates the initial lobe position
    void compute_initial_lobe_position( int idx_flow, Lobe & lobe );

    // perturbes the initial azimuthal angle of the lobe, which is
    void compute_lobe_axes( Lobe & lobe, const Vector2 & slope ) const; // computed from the terrain slope

    static void compute_descendent_lobe_position( Lobe & lobe, const Lobe & parent, Vector2 final_budding_point );

    void perturb_lobe_angle( Lobe & lobe, const Vector2 & slope );

    int select_parent_lobe( int idx_descendant );

    void compute_cumulative_descendents( std::vector<Lobe> & lobes ) const;

    void add_inertial_contribution( Lobe & lobe, const Lobe & parent, const Vector2 & slope ) const;

    void write_lobe_data_to_file( const std::vector<Lobe> & lobes, const std::filesystem::path & output_path );

    bool stop_condition( const Vector2 & point, double radius );

    void run();

private:
    int rng_seed;
    std::mt19937 gen{};
};

} // namespace Flowy