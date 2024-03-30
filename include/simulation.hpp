#pragma once
#include "asc_file.hpp"
#include "config.hpp"
#include "definitions.hpp"
#include "lobe.hpp"
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
    Simulation( const Config::InputParams & input )
            : input( input ),
              asc_file( AscFile( input.source ) ),
              lobe_dimensions( CommonLobeDimensions( input, asc_file ) ){};

    Config::InputParams input;
    AscFile asc_file;

    std::vector<Lobe> lobes; // Lobes per flows

    /*
Calculates the initial lobe position
*/
    void compute_initial_lobe_position( int idx_flow, int idx_lobe );

    void run();

private:
    int n_lobes_total = 0; // This is the total number of lobes, accumulated over all flows
    CommonLobeDimensions lobe_dimensions;
};

} // namespace Flowtastic