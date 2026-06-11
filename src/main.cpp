#include <omp.h>
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/config_parser.hpp"
#include "flowy/include/simulation.hpp"
#include "fmt/core.h"

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <argparse/argparse.hpp>
#include <filesystem>
namespace fs = std::filesystem;

int main( int argc, char * argv[] )
{
    using namespace Flowy;

    argparse::ArgumentParser program( "flowy" );

    program.add_argument( "config_file" ).help( "The config file to be used. Has to be in TOML format." );
    program.add_argument( "-a", "--asc_file" )
        .help( "The .asc file to be used for the terrain. This overwrites the `source` field in the input.toml file." );
    program.add_argument( "-n", "--name" )
        .help( "The run_name to be used. This overwrites the `run_name` in the input file and disables the system for "
               "automatically appending numbers to the run_name" );
    program.add_argument( "-o", "--output" )
        .help( fmt::format(
            "Specify the output directory. Defaults to `{}`", Config::InputParams().output_folder.string() ) );

    try
    {
        program.parse_args( argc, argv );
    }
    catch( const std::runtime_error & err )
    {
        fmt::print( stderr, "{}\n{}", err.what(), fmt::streamed( program ) );
        return 1;
    }

    fs::path config_file_path                      = program.get<std::string>( "config_file" );
    std::optional<fs::path> asc_file_path          = program.present<std::string>( "-a" );
    std::optional<std::string> output_dir_path_cli = program.present<std::string>( "-o" );
    std::optional<std::string> run_name            = program.present<std::string>( "-n" );

    auto input_params = Config::parse_config( config_file_path );
    validate_settings( input_params );

    if( asc_file_path.has_value() )
    {
        input_params.source = asc_file_path.value();
    }

    if( output_dir_path_cli.has_value() )
    {
        input_params.output_folder = output_dir_path_cli.value();
    }

    // lambda to get the name of the input backup file
    auto get_input_backup_name = [&]() { return fmt::format( "{}_inp.bak", input_params.run_name ); };

    // Decide run name
    // If the name is given through a CLI argument, we just accept it
    // and potentially overwrite files with the same run_name
    if( run_name.has_value() )
    {
        input_params.run_name = run_name.value();
    }
    else
    {
        // Else, the name comes from the input.toml file and we take it as a "base_name"
        // Then we check if the input backup file already exists for run_names with three appended digits:
        // That means we check for
        // {base_run_name}_000_inp.bak, {base_run_name}_001_inp.bak, ... etc.
        // As our final run name we take the first {base_run_name}_XXX that does not already have an
        // existing input backup file
        // (This is how Mr. Lava Loba handles the run_name)
        const std::string run_base_name  = input_params.run_name;
        constexpr int MAX_RUNS_IN_FOLDER = 999;
        for( int i = 0; i < MAX_RUNS_IN_FOLDER; i++ )
        {
            input_params.run_name = fmt::format( "{}_{:03}", run_base_name, i );
            if( !std::filesystem::is_regular_file( input_params.output_folder / get_input_backup_name() ) )
            {
                break;
            }
        }
    }

    // Back up input file
    // First we need to make sure that the output directory exists
    std::filesystem::create_directories( input_params.output_folder );
    auto path_to_inp_bak = input_params.output_folder / get_input_backup_name();

    // If the input backup already exists, we delete it
    if( std::filesystem::exists( path_to_inp_bak ) )
    {
        std::filesystem::remove( path_to_inp_bak );
    }
    // Then, we copy
    std::filesystem::copy_file( config_file_path, path_to_inp_bak );

    fmt::print( "=================================================================\n" );
    fmt::print( "Using input file: {}\n", config_file_path.string() );
    fmt::print( "Output directory path set to: {}\n", input_params.output_folder.string() );
    fmt::print( "run_name = {}\n", input_params.run_name );
    if( input_params.n_runs > 1 )
    {
        const int n_runs    = input_params.n_runs;
        const int base_seed = input_params.rng_seed.value_or( 0 );
        xt::xtensor<double, 2> sum;
        bool sum_init = false;
#pragma omp parallel
        {
            xt::xtensor<double, 2> local;
            bool local_init = false;
#pragma omp for schedule( dynamic )
            for( int i = 0; i < n_runs; i++ )
            {
                auto inp            = input_params;
                inp.write_lobes_csv = false;
                inp.output_folder   = input_params.output_folder / fmt::format( "ens_{}", i );
                inp.run_name        = "r";
                Simulation sim( inp, base_seed + i );
                sim.run();
                sim.compute_topography_thickness();
                if( !local_init )
                {
                    local      = sim.topography_thickness.height_data;
                    local_init = true;
                }
                else
                {
                    local += sim.topography_thickness.height_data;
                }
            }
#pragma omp critical
            {
                if( local_init )
                {
                    if( !sum_init )
                    {
                        sum      = local;
                        sum_init = true;
                    }
                    else
                    {
                        sum += local;
                    }
                }
            }
        }
        sum /= static_cast<double>( n_runs );
        auto grid = Simulation::construct_initial_topography( input_params );
        Flowy::Topography mean_topo( sum, grid.x_data, grid.y_data, DEFAULT_NO_DATA_VALUE_THICKNESS );
        auto out = input_params.output_folder / "ensemble_mean";
        std::filesystem::create_directories( input_params.output_folder );
        Flowy::AscFile( mean_topo, Flowy::OutputQuantity::Height ).save( out );
        fmt::print( "ensemble of {} runs -> {}\n", n_runs, out.string() );
    }
    else
    {
        auto simulation = Simulation( input_params, input_params.rng_seed );
        simulation.run();
    }
    fmt::print( "=================================================================\n" );
}
