#include "config_parser.hpp"
#include "fmt/core.h"
#include "simulation.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <argparse/argparse.hpp>
#include <filesystem>
namespace fs = std::filesystem;

int main( int argc, char * argv[] )
{
    using namespace Flowtastic;

    argparse::ArgumentParser program( "flowtastic" );

    program.add_argument( "config_file" ).help( "The config file to be used. Has to be in TOML format." );
    program.add_argument( "-a", "--asc_file" ).help( "The .asc file to be used for the terrain." );
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

    auto input_params = Config::parse_config( config_file_path );

    if( asc_file_path.has_value() )
    {
        input_params.source = asc_file_path.value();
    }

    if( output_dir_path_cli.has_value() )
    {
        input_params.output_folder = output_dir_path_cli.value();
    }

    fmt::print( "=================================================================\n" );
    fmt::print( "Using input file: {}\n", config_file_path.string() );
    fmt::print( "Output directory path set to: {}\n", input_params.output_folder.string() );
    auto simulation = Simulation( input_params, 0 );
    simulation.run();
    fmt::print( "=================================================================\n" );
}