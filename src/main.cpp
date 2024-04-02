#include "config_parser.hpp"
#include "simulation.hpp"

#include <fmt/format.h>
#include <fmt/ostream.h>
#include <argparse/argparse.hpp>
#include <filesystem>
namespace fs = std::filesystem;

int main( int argc, char * argv[] )
{
    argparse::ArgumentParser program( "flowtastic" );

    program.add_argument( "config_file" ).help( "The config file to be used. Has to be in TOML format." );

    try
    {
        program.parse_args( argc, argv );
    }
    catch( const std::runtime_error & err )
    {
        fmt::print( stderr, "{}\n{}", err.what(), fmt::streamed( program ) );
        return 1;
    }

    fs::path config_file_path = program.get<std::string>( "config_file" );

    auto input_params = Flowtastic::Config::parse_config( config_file_path );

    auto simulation   = Flowtastic::Simulation( input_params, 0 );

    simulation.run();
}