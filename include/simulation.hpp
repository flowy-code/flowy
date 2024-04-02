#pragma once
#include "asc_file.hpp"
#include "config.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
#include "topography.hpp"
#include <random>
#include <vector>

namespace Flowtastic
{

class CommonLobeDimensions
{
public:
    CommonLobeDimensions( const Config::InputParams & input, const AscFile & asc_file );

    double avg_lobe_thickness = 0;
    double lobe_area          = 0;
    double max_semiaxis       = 0;
    int max_cells             = 0;
    double thickness_min      = 0;
};

class Simulation
{
public:
    Simulation( const Config::InputParams & input, std::optional<int> rng_seed )
            : input( input ),
              asc_file( AscFile( input.source ) ),
              topography( Topography( asc_file ) ),
              lobe_dimensions( CommonLobeDimensions( input, asc_file ) ),
              gen( std::mt19937( rng_seed.value_or( std::random_device()() ) ) ){};

    Config::InputParams input;
    AscFile asc_file;
    Topography topography;

    std::vector<Lobe> lobes; // Lobes per flows

    // calculates the initial lobe position
    void compute_initial_lobe_position( int idx_flow, Lobe & lobe );

    // perturbes the initial azimuthal angle of the lobe, which is
    void compute_lobe_axes( Lobe & lobe, const Vector2 & slope ) const; // computed from the terrain slop
    void perturb_lobe_angle( Lobe & lobe, const Vector2 & slope );

    int select_parent_lobe( int idx_descendant );

    void run();

private:
    int n_lobes_total = 0; // This is the total number of lobes, accumulated over all flows
    CommonLobeDimensions lobe_dimensions;
    std::mt19937 gen{};
};

} // namespace Flowtastic