#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/asc_file.hpp"
#include "flowy/include/config.hpp"
#include "flowy/include/definitions.hpp"
#include "flowy/include/lobe.hpp"
#include "flowy/include/models/mr_lava_loba.hpp"
#include "flowy/include/topography.hpp"
#include "flowy/include/topography_file.hpp"
#include <chrono>
#include <filesystem>
#include <memory>
#include <optional>
#include <random>
#include <vector>

namespace Flowy
{

/**
 * @brief Track the current state of the Simulation run
 *
 */
struct SimulationState
{
    std::chrono::time_point<std::chrono::system_clock> t_run_start{};
    std::vector<Lobe> lobes{};

    int n_lobes_processed = 0;
    int n_lobes           = 0;
    int idx_flow          = 0;
    int idx_lobe          = 0;

    bool beginning_of_new_flow = true;
};

enum class FlowStatus
{
    Finished,
    Ongoing
};

enum class RunStatus
{
    Finished,
    Ongoing
};

class Simulation
{
public:
    Simulation( const Config::InputParams & input, std::optional<int> rng_seed );

    Config::InputParams input;

    Topography topography_initial;   // Stores the initial topography, before any simulations are run
    Topography topography_thickness; // Stores the height_difference between initial and final topography
    Topography topography;           // The topography, which is modified during the simulation
    std::vector<Lobe> lobes;         // Lobes per flows

    static Topography construct_initial_topography( const Config::InputParams & input );

    void compute_cumulative_descendents( std::vector<Lobe> & lobes ) const;

    void write_lobe_data_to_file( const std::vector<Lobe> & lobes, const std::filesystem::path & output_path );

    bool stop_condition( const Vector2 & point, double radius ) const;

    void write_avg_thickness_file();

    // Check if the dem has to be written (because of the input.write_dem_every_n_lobes_setting) and, if yes, writes the
    // topography
    void write_thickness_if_necessary( int n_lobes_processed );

    // Computes the topography_thickness field by subtracting the initial topography and dividing by (1.0 - filling_parameter)
    void compute_topography_thickness();

    std::unique_ptr<TopographyFile>
    get_file_handle( const Topography & topography, OutputQuantity output_quantity ) const;

    void save_post_run_output();

    void run();

    // Perform `n_steps` steps of the simulation (per step, a single lobe is added to the topography)
    RunStatus steps( int n_steps );

    void reset_simulation_state()
    {
        simulation_state = std::nullopt;
    }

    std::optional<SimulationState> get_simulation_state() const
    {
        return simulation_state;
    }

private:
    int rng_seed;
    std::mt19937 gen{};
    std::optional<SimulationState> simulation_state = std::nullopt;

    void post_flow_hook( int idx_flow, std::vector<Lobe> & lobes );
    void process_initial_lobe( int idx_flow, Lobe & lobe_cur, const CommonLobeDimensions & common_lobe_dimensions );
    FlowStatus process_descendent_lobe(
        int idx_lobe, std::vector<Lobe> & lobes, const CommonLobeDimensions & common_lobe_dimensions );
};

} // namespace Flowy
