#include "config_parser.hpp"
#include "config.hpp"
#include <toml++/toml.h>
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

    // our own
    set_if_specified( params.write_lobes_csv, tbl["write_lobes_csv"] );
    set_if_specified( params.print_remaining_time, tbl["print_remaining_time"] );
    params.rng_seed = tbl["rng_seed"].value<int>();

    // From input.py
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

    params.fissure_probabilities = tbl["fissure_probabilities"].value<double>();
    params.total_volume          = tbl["total_volume"].value<double>();
    params.east_to_vent          = tbl["east_to_vent"].value<double>();
    params.west_to_vent          = tbl["west_to_vent"].value<double>();
    params.south_to_vent         = tbl["south_to_vent"].value<double>();
    params.north_to_vent         = tbl["north_to_vent"].value<double>();
    params.channel_file          = tbl["channel_file"].value<std::string>();
    params.alfa_channel          = tbl["alfa_channel"].value<double>();
    params.d1                    = tbl["d1"].value<double>();
    params.d2                    = tbl["d2"].value<double>();
    params.eps                   = tbl["eps"].value<double>();
    params.union_diff_file       = tbl["union_diff_file"].value<std::string>();

    // from input_advanced.py
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

    auto restart_files                = parse_vector<std::string>( tbl["Advanced"]["restart_files"] );
    params.restart_filling_parameters = parse_vector<double>( tbl["Advanced"]["restart_filling_parameters"] );

    if( !restart_files.empty() )
    {
        params.restart_files = std::vector<std::filesystem::path>{};
    }
    for( auto & s : restart_files )
    {
        params.restart_files.value().emplace_back( s );
    }

    return params;
}
} // namespace Flowy::Config
