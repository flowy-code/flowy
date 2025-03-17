// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/simulation.hpp"
#include "flowy/include/asc_file.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/math.hpp"
#include "flowy/include/models/mr_lava_loba.hpp"
#include "flowy/include/topography.hpp"
#include "flowy/include/topography_file.hpp"
#include "fmt/core.h"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>

// This is more for clarity since the ifdef is in the netcdf_file as well
#ifdef WITH_NETCDF
#include "flowy/include/netcdf_file.hpp"
#endif

namespace Flowy
{

Simulation::Simulation( const Config::InputParams & input, std::optional<int> rng_seed ) : input( input )
{
    this->rng_seed = rng_seed.value_or( std::random_device()() );
    gen            = std::mt19937( this->rng_seed );

    // Create output directory
    std::filesystem::create_directories( input.output_folder ); // Create the output directory

    topography = Simulation::construct_initial_topography( input );

    // Make a copy of the initial topography
    topography_initial = topography;
}

Topography Simulation::construct_initial_topography( const Config::InputParams & input )
{
    AscFile asc_file{};

    // Crop if all of these have a value
    if( input.east_to_vent.has_value() && input.west_to_vent.has_value() && input.south_to_vent.has_value()
        && input.north_to_vent.has_value() )
    {
        auto crop = TopographyCrop{};

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

    auto topography = asc_file.to_topography();

    if( input.restart_files.has_value() )
    {
        auto n_restart = input.restart_files->size();

        if( input.restart_filling_parameters.has_value() && input.restart_filling_parameters->size() != n_restart )
        {
            throw std::runtime_error( fmt::format(
                "restart_filling_parameters has size {}, but there are {} restart files",
                input.restart_filling_parameters->size(), n_restart ) );
        }

        for( size_t i_restart = 0; i_restart < n_restart; i_restart++ )
        {
            double filling_parameter = 1.0;
            if( input.restart_filling_parameters.has_value() )
            {
                filling_parameter = input.restart_filling_parameters.value()[i_restart];
            }
            fmt::print(
                "Restart file {}/{}: {}\n", i_restart + 1, n_restart, input.restart_files.value()[i_restart].string() );
            fmt::print( "    filling_parameter = {}\n", filling_parameter );

            // Create the topography for the restart file
            auto asc_file_restart   = AscFile( input.restart_files.value()[i_restart] );
            auto restart_topography = asc_file_restart.to_topography();
            // Add the restart file to the topography
            topography.add_to_topography( restart_topography, filling_parameter );
        }
    }

    return topography;
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

// Depth first search to compute cumulative descendents
// TODO for flows with a very large number of lobes, the recursion might become a problem
// It would be better to write a recrusion free version of this function.
int dfs( int lobe_idx, std::vector<Lobe> & lobes, std::vector<std::vector<int>> & child_node_list )
{
    int total_descendants = 0;
    auto & child_nodes    = child_node_list[lobe_idx];
    for( int child : child_nodes )
    {
        total_descendants += 1 + dfs( child, lobes, child_node_list );
    }
    lobes[lobe_idx].n_descendents = total_descendants;
    return total_descendants;
}

void Simulation::compute_cumulative_descendents( std::vector<Lobe> & lobes ) const
{
    // First we invert the parent-child relationship by recording a list of child node indices for each lobe
    std::vector<std::vector<int>> child_node_list( lobes.size() );
    for( size_t i_lobe = 0; i_lobe < lobes.size(); i_lobe++ )
    {
        const Lobe & cur_lobe = lobes[i_lobe];
        if( cur_lobe.idx_parent.has_value() )
        {
            child_node_list[cur_lobe.idx_parent.value()].push_back( i_lobe );
        }
    }

    // Then, we have to start a depth first search separately on each root
    for( int i_root = 0; i_root < input.n_init; i_root++ )
    {
        dfs( i_root, lobes, child_node_list );
    }
}

bool Simulation::stop_condition( const Vector2 & point, double radius ) const
{
    return topography.is_point_near_boundary( point, radius )
           || topography.get_height( point ) <= topography.no_data_value;
}

void Simulation::write_avg_thickness_file()
{
    const auto path = input.output_folder / fmt::format( "{}_avg_thick.txt", input.run_name );

    std::fstream file;
    file.open( path, std::fstream::in | std::fstream::out | std::fstream::trunc );

    if( !file.is_open() )
    {
        throw std::runtime_error( fmt::format( "Unable to create file: '{}'", path.string() ) );
    }

    const double volume        = topography_thickness.volume();
    const double area          = topography_thickness.area();
    const double avg_thickness = volume / area;
    const double cell_area     = topography.cell_size() * topography.cell_size();

    auto lobe_dimensions = CommonLobeDimensions( input );

    file << fmt::format( "Average lobe thickness = {} m\n", lobe_dimensions.avg_lobe_thickness );
    file << fmt::format( "Total volume = {} m3\n", volume );
    file << fmt::format( "Total area = {} m2\n", area );
    file << fmt::format( "Average thickness full = {} m\n", avg_thickness );

    // Create a flattened, sorted view of the thickness, which will be used in the bisection search later
    auto thickness_non_zero = xt::filter( topography_thickness.height_data, topography_thickness.height_data > 0 );

    if( thickness_non_zero.size() == 0 )
    {
        fmt::print( "Cannot determine masked quantities, since thickness_non_zero.size() == 0\n" );
        file.close();
        return;
    }

    auto flatten          = xt::flatten( thickness_non_zero );
    auto thickness_sorted = xt::eval( xt::sort( flatten ) );
    const int n_cells     = thickness_sorted.size();

    // This lambda performs bisection search to find the threshold thickness at which a
    // relative volume proportion of `thresh` is contained within cells with greater thickness than the threshold thickness
    auto bisection_search = [&]( double thresh, double tol, int max_iter )
    {
        int idx_lo = 0;
        int idx_hi = n_cells - 1;

        // If the relation between masked volume and threshold thickness would be linear,
        // this would be the solution for the index. Therefore, we use it as an initial guess.
        int idx_cur = std::max<int>( ( n_cells - 1 ) * ( 1.0 - thresh ), 1 );

        double volume_cur{};
        double ratio{};
        for( int iter = 0; iter < max_iter; iter++ )
        {
            volume_cur = cell_area * xt::sum( xt::view( thickness_sorted, xt::range( idx_cur, -1 ) ) )();

            // The ratio between the sum of the flow values is the same as the volume ratio,
            // since the cell_size cancels out
            ratio = volume_cur / volume;

            // Stop if we are within tol
            if( std::abs( volume_cur / volume - thresh ) < tol )
            {
                break;
            }

            if( ratio > thresh )
            {
                idx_lo = idx_cur;
            }
            else
            {
                idx_hi = idx_cur;
            }

            idx_cur = std::clamp<int>( 0.5 * ( idx_lo + idx_hi ), 0, thickness_sorted.size() - 1 );

            // fmt::print( "iter {} idx {} ratio {}\n", iter, idx_cur, ratio );
        }

        double threshold_thickness = thickness_sorted[idx_cur];
        double area_cur = cell_area * xt::count_nonzero( xt::view( thickness_sorted, xt::range( idx_cur, -1 ) ) )();

        return std::array<double, 4>{ threshold_thickness, volume_cur, area_cur, ratio };
    };

    for( auto & threshold : input.masking_threshold )
    {
        auto const [threshold_thickness, volume_masked, area_masked, ratio]
            = bisection_search( threshold, input.masking_tolerance, input.masking_max_iter );

        double avg_thickness = volume_masked / area_masked;

        file << fmt::format( "Masking threshold = {}\n", threshold );
        file << fmt::format( "Masked volume = {} m3\n", volume_masked );
        file << fmt::format( "Masked area = {} m2\n", area_masked );
        file << fmt::format( "Average thickness mask = {} m\n", avg_thickness );

        // We make a temporary copy of the thickness topography
        Topography topography_masked = topography_thickness;

        // Apply the thickness filter to both
        xt::filter( topography_masked.height_data, topography_thickness.height_data < threshold_thickness ) = 0.0;
        xt::filter( topography_masked.hazard, topography_thickness.height_data < threshold_thickness )      = 0.0;

        // Write the masked thickness and the masked hazard maps
        auto file_thick = get_file_handle( topography_masked, OutputQuantity::Height );
        file_thick->save(
            input.output_folder / fmt::format( "{}_thickness_masked_{:.2f}", input.run_name, threshold ) );

        if( input.save_hazard_data )
        {
            auto file_hazard = get_file_handle( topography_masked, OutputQuantity::Hazard );
            file_hazard->save(
                input.output_folder / fmt::format( "{}_hazard_masked_{:.2f}", input.run_name, threshold ) );
        }
    }
    file.close();
}

std::unique_ptr<TopographyFile>
Simulation::get_file_handle( const Topography & topography, OutputQuantity output_quantity ) const
{
    std::unique_ptr<TopographyFile> res{};

    if( input.output_settings.use_netcdf )
    {
#ifdef WITH_NETCDF
        auto netcdf_file        = NetCDFFile( topography, output_quantity );
        netcdf_file.compression = input.output_settings.compression;
        netcdf_file.data_type   = input.output_settings.data_type;
        if( netcdf_file.compression )
        {
            netcdf_file.compression_level = input.output_settings.compression_level;
            netcdf_file.shuffle           = input.output_settings.shuffle;
        }
        res = std::make_unique<NetCDFFile>( netcdf_file );
#else
        throw std::runtime_error( "NetCDF requested but Flowy wasn't built with NetCDF support!\n" );
#endif
    }
    else
    {
        auto asc_file = AscFile( topography, output_quantity );
        res           = std::make_unique<AscFile>( asc_file );
    }

    if( input.output_settings.crop_to_content )
        res->crop_to_content();

    return res;
}

void Simulation::run()
{
    // Initialize MrLavaLoba method
    auto mr_lava_loba = MrLavaLoba( input, gen );

    int n_lobes_processed = 0;

    // Make a copy of the initial topography
    auto t_run_start = std::chrono::high_resolution_clock::now();

    for( int idx_flow = 0; idx_flow < input.n_flows; idx_flow++ )
    {
        // Determine n_lobes
        int n_lobes = mr_lava_loba.compute_n_lobes( idx_flow );

        lobes = std::vector<Lobe>{};
        lobes.reserve( n_lobes );

        // set the intersection cache
        topography.reset_intersection_cache( n_lobes );

        // Build initial lobes which do not propagate descendents
        for( int idx_lobe = 0; idx_lobe < input.n_init; idx_lobe++ )
        {
            lobes.emplace_back();
            Lobe & lobe_cur = lobes.back();

            mr_lava_loba.compute_initial_lobe_position( idx_flow, lobe_cur );

            // Compute the thickness of the lobe
            lobe_cur.thickness = mr_lava_loba.compute_current_lobe_thickness( idx_lobe, n_lobes );

            auto [height_lobe_center, slope] = topography.height_and_slope( lobe_cur.center );

            if( height_lobe_center == topography.no_data_value )
            {
                throw std::runtime_error(
                    "The initial lobe center has been placed on a no_data value point in the topography." );
            }

            // Perturb the angle (and set it)
            lobe_cur.set_azimuthal_angle( std::atan2( slope[1], slope[0] ) ); // Sets the angle prior to perturbation
            const double slope_norm = xt::linalg::norm( slope, 2 );           // Similar to np.linalg.norm
            mr_lava_loba.perturb_lobe_angle( lobe_cur, slope_norm );

            // compute lobe axes
            mr_lava_loba.compute_lobe_axes( lobe_cur, slope_norm );

            // Add rasterized lobe
            topography.add_lobe( lobe_cur, input.volume_correction, idx_lobe );
            n_lobes_processed++;
        }

        // Loop over the rest of the lobes (skipping the initial ones).
        // Each lobe is a descendant of a parent lobe
        for( int idx_lobe = input.n_init; idx_lobe < n_lobes; idx_lobe++ )
        {
            lobes.emplace_back();
            Lobe & lobe_cur = lobes.back();

            // Select which of the previously created lobes is the parent lobe
            // from which the new descendent lobe will bud
            auto idx_parent    = mr_lava_loba.select_parent_lobe( idx_lobe, lobes );
            Lobe & lobe_parent = lobes[idx_parent];

            // stopping condition (parent lobe close the domain boundary or at a not defined z value)
            if( stop_condition( lobe_parent.center, lobe_parent.semi_axes[0] ) )
            {
                lobes.pop_back();
                break;
            }

            // Find the preliminary budding point on the perimeter of the parent lobe (npoints is the number of raster
            // points on the ellipse)
            Flowy::Vector2 budding_point = topography.find_preliminary_budding_point( lobe_parent, input.npoints );

            auto [height_lobe_center, slope_parent] = topography.height_and_slope( lobe_parent.center );
            auto [height_bp, slope_bp]              = topography.height_and_slope( budding_point );

            Vector2 diff = ( budding_point - lobe_parent.center );

            // Perturb the angle and set it (not on the parent anymore)
            lobe_cur.set_azimuthal_angle( std::atan2( diff[1], diff[0] ) ); // Sets the angle prior to perturbation
            const double slope_parent_norm = topography.slope_between_points( lobe_parent.center, budding_point );
            mr_lava_loba.perturb_lobe_angle( lobe_cur, slope_parent_norm );

            // Add the inertial contribution
            mr_lava_loba.add_inertial_contribution( lobe_cur, lobe_parent, slope_parent_norm );

            // Compute the final budding point
            // It is defined by the point on the perimeter of the parent lobe closest to the center of the new lobe
            auto angle_diff             = lobe_parent.get_azimuthal_angle() - lobe_cur.get_azimuthal_angle();
            Vector2 final_budding_point = lobe_parent.point_at_angle( -angle_diff );

            // final_budding_point = budding_point;
            if( stop_condition( final_budding_point, lobe_parent.semi_axes[0] ) )
            {
                lobes.pop_back();
                break;
            }
            // Get the slope at the final budding point
            double slope_budding_point = topography.slope_between_points( lobe_parent.center, final_budding_point );

            // compute the new lobe axes
            mr_lava_loba.compute_lobe_axes( lobe_cur, slope_budding_point );

            // Get new lobe center
            mr_lava_loba.compute_descendent_lobe_position( lobe_cur, lobe_parent, final_budding_point );

            if( stop_condition( lobe_cur.center, lobe_cur.semi_axes[0] ) )
            {
                lobes.pop_back();
                break;
            }

            // Compute the thickness of the lobe
            lobe_cur.thickness = mr_lava_loba.compute_current_lobe_thickness( idx_lobe, n_lobes );

            // Add rasterized lobe
            topography.add_lobe( lobe_cur, input.volume_correction, idx_lobe );
            n_lobes_processed++;
        }

        if( input.save_hazard_data )
        {
            compute_cumulative_descendents( lobes );
            topography.compute_hazard_flow( lobes );
        }

        if( input.write_lobes_csv )
        {
            write_lobe_data_to_file( lobes, input.output_folder / fmt::format( "lobes_{}.csv", idx_flow ) );
        }

        if( input.print_remaining_time )
        {
            auto t_cur          = std::chrono::high_resolution_clock::now();
            auto remaining_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                ( input.n_flows - idx_flow - 1 ) * ( t_cur - t_run_start ) / ( idx_flow + 1 ) );
            fmt::print( "     remaining_time = {:%Hh %Mm %Ss}\n", remaining_time );
        }
    }

    auto t_cur      = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>( ( t_cur - t_run_start ) );
    fmt::print( "total_time = {:%Hh %Mm %Ss}\n", total_time );

    fmt::print( "Total number of processed lobes = {}\n", n_lobes_processed );

    if( total_time.count() > 0 )
    {
        auto lobes_per_ms = n_lobes_processed / total_time.count();
        fmt::print( "n_lobes/ms = {}\n", lobes_per_ms );
    }

    fmt::print( "Used RNG seed: {}\n", rng_seed );

    // Save initial topography to asc file
    auto file_initial = get_file_handle( topography_initial, OutputQuantity::Height );
    file_initial->save( input.output_folder / fmt::format( "{}_DEM", input.run_name ) );

    // Save final topography to asc file
    if( input.save_final_dem )
    {
        auto file_final = get_file_handle( topography, OutputQuantity::Height );
        file_final->save( input.output_folder / fmt::format( "{}_DEM_final", input.run_name ) );
    }

    // Save full thickness to asc file
    topography_thickness               = topography;
    topography_thickness.no_data_value = DEFAULT_NO_DATA_VALUE_THICKNESS;
    topography_thickness.height_data -= topography_initial.height_data;
    topography_thickness.height_data /= ( 1.0 - input.thickening_parameter );

    auto file_thick = get_file_handle( topography_thickness, OutputQuantity::Height );
    file_thick->save( input.output_folder / fmt::format( "{}_thickness_full", input.run_name ) );

    // Save the full hazard map
    if( input.save_hazard_data )
    {
        auto file_hazard = get_file_handle( topography, OutputQuantity::Hazard );
        file_hazard->save( input.output_folder / fmt::format( "{}_hazard_full", input.run_name ) );
    }

    write_avg_thickness_file();
}

} // namespace Flowy
