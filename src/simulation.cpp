#include "simulation.hpp"
#include "math.hpp"
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

void Simulation::compute_initial_lobe_position( int idx_flow, int idx_lobe )
{
    // For now, we've just implemented vent_flag = 0
    int idx_vent           = std::floor( idx_flow * input.n_vents() / input.n_flows );
    lobes[idx_lobe].center = input.vent_coordinates[idx_vent];
}

void Simulation::run()
{
    for( int idx_flow = 0; idx_flow < input.n_flows; idx_flow++ )
    {
        // Determine n_lobes
        // @TODO: Calculate this later, using a beta distribution etc.
        int n_lobes = 0.5 * ( input.min_n_lobes + input.max_n_lobes );
        n_lobes_total += n_lobes;

        double delta_lobe_thickness
            = 2.0 * ( lobe_dimensions.avg_lobe_thickness - lobe_dimensions.thickness_min ) / ( n_lobes - 1.0 );

        // Build initial lobes which do  not propagate descendents
        for( int idx_lobe = 0; idx_lobe < input.n_init; idx_lobe++ )
        {
            compute_initial_lobe_position( idx_lobe, idx_flow );
            // get the slope
            // perturb the angle
            // compute lobe axes
            // rasterize_lobe
        }
    }
}

} // namespace Flowtastic