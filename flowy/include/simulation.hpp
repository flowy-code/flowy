#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/asc_file.hpp"
#include "flowy/include/config.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/topography.hpp"
#include <filesystem>
#include <random>
#include <vector>

namespace Flowy
{

class CommonLobeDimensions
{
public:
    CommonLobeDimensions() = default;
    CommonLobeDimensions( const Config::InputParams & input, const AscFile & asc_file );

    double avg_lobe_thickness = 0;
    double lobe_area          = 0;
    double max_semiaxis       = 0;
    int max_cells             = 0;
    double thickness_min      = 0;
    double exp_lobe_exponent  = 1;
};

class Simulation
{
public:
    Simulation( const Config::InputParams & input, std::optional<int> rng_seed );

    Config::InputParams input;
    AscFile asc_file;
    Topography topography_initial;   // Stores the initial topography, before any simulations are run
    Topography topography_thickness; // Stores the height_difference between initial and final topography
    Topography topography;           // The topography, which is modified during the simulation
    CommonLobeDimensions lobe_dimensions;

    std::vector<Lobe> lobes; // Lobes per flows

    // calculates the initial lobe position
    void compute_initial_lobe_position( int idx_flow, Lobe & lobe );

    // perturbes the initial azimuthal angle of the lobe, which is
    void compute_lobe_axes( Lobe & lobe, const Vector2 & slope ) const; // computed from the terrain slope

    void compute_descendent_lobe_position( Lobe & lobe, const Lobe & parent, Vector2 final_budding_point );

    void perturb_lobe_angle( Lobe & lobe, double slope );

    int select_parent_lobe( int idx_descendant );

    void compute_cumulative_descendents( std::vector<Lobe> & lobes ) const;

    void add_inertial_contribution( Lobe & lobe, const Lobe & parent, double slope ) const;

    void write_lobe_data_to_file( const std::vector<Lobe> & lobes, const std::filesystem::path & output_path );

    bool stop_condition( const Vector2 & point, double radius );

    void write_avg_thickness_file();

    std::optional<std::vector<double>> compute_cumulative_fissure_length();

    void run();

private:
    int rng_seed;
    std::mt19937 gen{};
};

} // namespace Flowy
