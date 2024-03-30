#pragma once
#include "asc_file.hpp"
#include "config.hpp"

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

private:
    CommonLobeDimensions lobe_dimensions;
};

} // namespace Flowtastic