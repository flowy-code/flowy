// GPL v3 License
// Copyright 2023--present Flowy developers

#include "flowy/include/config_parser.hpp"
#include "flowy/include/config.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/netcdf_file.hpp"
#include "thirdparty/toml.hpp"
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>
#include <filesystem>
#include <optional>
#include <stdexcept>

namespace Flowy::Config
{

void set_if_specified( auto & opt, const auto & toml_opt )
{
    using T    = typename std::remove_reference<decltype( opt )>::type;
    auto t_opt = toml_opt.template value<T>();
    if( t_opt.has_value() )
        opt = t_opt.value();
}

template<typename T>
std::vector<T> parse_vector( const auto & toml_node )
{
    std::vector<T> res{};

    if( toml_node.is_array() )
    {
        toml::array * toml_arr = toml_node.as_array();
        toml_arr->for_each( [&]( auto && elem ) { res.push_back( elem.value_or( T{} ) ); } );
    }

    return res;
}

InputParams parse_config( const std::filesystem::path & path )
{
    InputParams params;

    toml::table tbl;
    tbl = toml::parse_file( path.string() );

    // ===================================================================================
    // new settings
    // ===================================================================================

    set_if_specified( params.write_lobes_csv, tbl["write_lobes_csv"] );
    set_if_specified( params.print_remaining_time, tbl["print_remaining_time"] );
    set_if_specified( params.save_final_dem, tbl["save_final_dem"] );

    std::optional<std::string> output_folder_string{};
    output_folder_string = tbl["output_folder"].value<std::string>();
    if( output_folder_string.has_value() )
    {
        params.output_folder = output_folder_string.value();
    }

    params.rng_seed = tbl["rng_seed"].value<int>();

    set_if_specified( params.masking_tolerance, tbl["masking_tolerance"] );
    set_if_specified( params.masking_max_iter, tbl["masking_max_iter"] );
    set_if_specified( params.volume_correction, tbl["volume_correction"] );

    // Output
    set_if_specified( params.output_settings.crop_to_content, tbl["Output"]["crop_to_content"] );
    set_if_specified( params.output_settings.use_netcdf, tbl["Output"]["use_netcdf"] );
    set_if_specified( params.output_settings.compression, tbl["Output"]["compression"] );
    set_if_specified( params.output_settings.compression_level, tbl["Output"]["compression_level"] );
    set_if_specified( params.output_settings.shuffle, tbl["Output"]["shuffle"] );

    std::optional<std::string> packing_data_type;
    packing_data_type = tbl["Output"]["packing_data_type"].value<std::string>();

    if( packing_data_type == "float" )
    {
        params.output_settings.data_type = StorageDataType::Float;
    }
    else if( packing_data_type == "double" )
    {
        params.output_settings.data_type = StorageDataType::Double;
    }
    else if( packing_data_type == "short" )
    {
        params.output_settings.data_type = StorageDataType::Short;
    }
    else if( packing_data_type.has_value() )
    {
        throw std::runtime_error( fmt::format( "Unknown packing_data_type: '{}'", packing_data_type.value() ) );
    }

    // set_if_specified( params.output_settings.data_type, tbl["Output"]["Packing_data_type"] );

    // ===================================================================================
    // mr lava loba settings from input.py
    // ===================================================================================
    set_if_specified( params.run_name, tbl["run_name"] );

    std::string source_string;
    set_if_specified( source_string, tbl["source"] );
    params.source = std::filesystem::path( source_string );

    std::vector<double> x_vent = parse_vector<double>( tbl["x_vent"] );
    std::vector<double> y_vent = parse_vector<double>( tbl["y_vent"] );
    if( x_vent.size() != y_vent.size() )
    {
        throw std::runtime_error( "x_vent and y_vent have different sizes" );
    }
    for( size_t i = 0; i < x_vent.size(); i++ )
    {
        params.vent_coordinates.push_back( { x_vent[i], y_vent[i] } );
    }

    std::vector<double> x_vent_end = parse_vector<double>( tbl["x_vent_end"] );
    std::vector<double> y_vent_end = parse_vector<double>( tbl["y_vent_end"] );

    if( x_vent_end.size() != y_vent_end.size() )
    {
        throw std::runtime_error( "x_vent_end and y_vent_end have different sizes" );
    }

    if( !x_vent_end.empty() )
    {
        std::vector<Vector2> fissure_end_coordinates{};
        for( size_t i = 0; i < x_vent_end.size(); i++ )
        {
            fissure_end_coordinates.push_back( { x_vent_end[i], y_vent_end[i] } );
        }
        params.fissure_end_coordinates = fissure_end_coordinates;
    }

    std::optional<int> hazard_flag = tbl["hazard_flag"].value<int>();
    if( hazard_flag.has_value() )
    {
        params.save_hazard_data = ( hazard_flag.value() == 1 );
    }

    // Masking threshold can be given as a single float like
    //      masking_threshold = 0.97
    // or as an array like
    //      masking_threshold = [0.97, 0.95, 0.92]
    // therefore, we need special logic to parse it
    if( tbl["masking_threshold"].is_array() )
    {
        params.masking_threshold = parse_vector<double>( tbl["masking_threshold"] );
    }
    else
    {
        std::optional<double> val = tbl["masking_threshold"].value<double>();
        if( val.has_value() )
        {
            params.masking_threshold.push_back( val.value() );
        }
    }

    set_if_specified( params.n_flows, tbl["n_flows"] );
    set_if_specified( params.n_lobes, tbl["n_lobes"] );
    set_if_specified( params.thickening_parameter, tbl["thickening_parameter"] );
    params.prescribed_lobe_area          = tbl["lobe_area"].value<double>();
    params.prescribed_avg_lobe_thickness = tbl["avg_lobe_thickness"].value<double>();

    set_if_specified( params.min_n_lobes, tbl["min_n_lobes"] );
    set_if_specified( params.max_n_lobes, tbl["max_n_lobes"] );
    set_if_specified( params.inertial_exponent, tbl["inertial_exponent"] );
    set_if_specified( params.lobe_exponent, tbl["lobe_exponent"] );
    set_if_specified( params.max_slope_prob, tbl["max_slope_prob"] );
    set_if_specified( params.thickness_ratio, tbl["thickness_ratio"] );
    set_if_specified( params.fixed_dimension_flag, tbl["fixed_dimension_flag"] );

    set_if_specified( params.vent_flag, tbl["vent_flag"] );

    auto fissure_probabilities = parse_vector<double>( tbl["fissure_probabilities"] );
    if( !fissure_probabilities.empty() )
        params.fissure_probabilities = fissure_probabilities;

    params.total_volume    = tbl["total_volume"].value<double>();
    params.east_to_vent    = tbl["east_to_vent"].value<double>();
    params.west_to_vent    = tbl["west_to_vent"].value<double>();
    params.south_to_vent   = tbl["south_to_vent"].value<double>();
    params.north_to_vent   = tbl["north_to_vent"].value<double>();
    params.channel_file    = tbl["channel_file"].value<std::string>();
    params.alfa_channel    = tbl["alfa_channel"].value<double>();
    params.d1              = tbl["d1"].value<double>();
    params.d2              = tbl["d2"].value<double>();
    params.eps             = tbl["eps"].value<double>();
    params.union_diff_file = tbl["union_diff_file"].value<std::string>();

    // ===================================================================================
    // mr lava loba settings from input_advanced.py
    // ===================================================================================
    set_if_specified( params.npoints, tbl["Advanced"]["npoints"] );
    set_if_specified( params.n_init, tbl["Advanced"]["n_init"] );
    set_if_specified( params.dist_fact, tbl["Advanced"]["dist_fact"] );
    set_if_specified( params.flag_threshold, tbl["Advanced"]["flag_threshold"] );
    set_if_specified( params.a_beta, tbl["Advanced"]["a_beta"] );
    set_if_specified( params.b_beta, tbl["Advanced"]["b_beta"] );
    set_if_specified( params.max_aspect_ratio, tbl["Advanced"]["max_aspect_ratio"] );
    set_if_specified( params.saveraster_flag, tbl["Advanced"]["saveraster_flag"] );
    set_if_specified( params.aspect_ratio_coeff, tbl["Advanced"]["aspect_ratio_coeff"] );
    set_if_specified( params.start_from_dist_flag, tbl["Advanced"]["start_from_dist_flag"] );
    set_if_specified( params.force_max_length, tbl["Advanced"]["force_max_length"] );
    set_if_specified( params.max_length, tbl["Advanced"]["max_length"] );
    set_if_specified( params.n_check_loop, tbl["Advanced"]["n_check_loop"] );

    auto restart_files = parse_vector<std::string>( tbl["Advanced"]["restart_files"] );
    if( !restart_files.empty() )
    {
        params.restart_files = std::vector<std::filesystem::path>{};
    }
    for( auto & s : restart_files )
    {
        params.restart_files.value().emplace_back( s );
    }

    auto restart_filling_parameters = parse_vector<double>( tbl["Advanced"]["restart_filling_parameters"] );
    if( !restart_filling_parameters.empty() )
    {
        params.restart_filling_parameters = restart_filling_parameters;
    }

    return params;
}

// This macro expands the variable x to: "x", x and is used together with the check function
#define name_and_var( x ) #x, x

// Helper function to check variables, depending on a condition function with an optional explanation
void check(
    const std::string & variable_name, auto variable, auto condition,
    const std::optional<std::string> & explanation = std::nullopt )
{
    if( !condition( variable ) )
    {
        std::string msg = fmt::format( "The value {} is not valid for {}", variable, variable_name );
        if( explanation.has_value() )
        {
            msg += "\n";
            msg += explanation.value();
        }
        throw std::runtime_error( msg );
    }
}

void validate_settings( const InputParams & options )
{
    auto g_zero           = []( auto x ) { return x > 0; };
    auto geq_zero         = []( auto x ) { return x >= 0; };
    auto geq_zero_leq_one = []( auto x ) { return x >= 0 && x <= 1; };

    check( name_and_var( options.n_flows ), g_zero );
    check( name_and_var( options.min_n_lobes ), geq_zero );
    check( name_and_var( options.max_n_lobes ), g_zero );
    check( name_and_var( options.lobe_exponent ), geq_zero_leq_one );
    check( name_and_var( options.max_slope_prob ), geq_zero_leq_one );
    check( name_and_var( options.inertial_exponent ), geq_zero );
    check( name_and_var( options.a_beta ), geq_zero );
    check( name_and_var( options.b_beta ), geq_zero );
    check( name_and_var( options.n_init ), []( auto x ) { return x >= 1; } );
    check( name_and_var( options.dist_fact ), geq_zero );
    check( name_and_var( options.npoints ), []( auto x ) { return x >= 1; } );
    check( name_and_var( options.aspect_ratio_coeff ), geq_zero );
    check( name_and_var( options.max_aspect_ratio ), g_zero );
    check( name_and_var( options.thickening_parameter ), []( auto x ) { return x >= 0.0 && x < 1.0; } );
    check(
        name_and_var( options.vent_flag ), []( auto x ) { return x >= 0.0 && x < 9; },
        "Allowed values of the vent_flag are 0 to 8, inclusive." );

    // Output settings validation
    check(
        name_and_var( options.output_settings.compression_level ), []( auto x ) { return x >= 0 && x <= 9; },
        "The compression level can only be between 0 and 9, inclusive." );

    check(
        name_and_var( options.fissure_end_coordinates ),
        [&]( auto & x )
        {
            if( x.has_value() )
            {
                return options.vent_coordinates.size() == x->size();
            }
            return true;
        },
        "x_vent and x_vent_end have different sizes" );
}

} // namespace Flowy::Config
