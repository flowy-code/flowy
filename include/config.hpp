#pragma once 
#include <string>
#include <filesystem>
#include "definitions.hpp"

namespace Flowtastic::Config {

    class InputParams{
        public:
    std::filesystem::path run_name{}; // Name of the run (used to save the parameters and the output)
    std::filesystem::path source{}; // File name of ASCII digital elevation model (.asc file)
    MatrixX vent_coordinates{}; // of shape [n_vents, 2]
    bool save_hazard_data{}; // If true, a raster map is saved, such that the values represent the probability of a cell to be covered by lava.
    int n_flows{}; // Number of flows per run
    int n_lobes{}; // Number of lobes per flow

    // Variables we don't understand
    // masking_threshold: float = 0
    // min_n_lobes: int = 0
    // max_n_lobes: int = 0
    thickening_parameter: float = 0
    lobe_area: float = 0
    inertial_exponent: float = 0
    lobe_exponent: float = 0
    max_slope_prob: float = 0
    thickness_ratio: float = 0
    fixed_dimension_flag: int = 0
    vent_flag: int = 0
    fissure_probabilities: Optional[float] = None
    total_volume: Optional[float] = None
    volume_flag: Optional[int] = None
    east_to_vent: Optional[float] = None
    west_to_vent: Optional[float] = None
    south_to_vent: Optional[float] = None
    north_to_vent: Optional[float] = None
    channel_file: Optional[str] = None
    alfa_channel: Optional[float] = None
    d1: Optional[float] = None
    d2: Optional[float] = None
    eps: Optional[float] = None
    union_diff_file: Optional[str] = None

    # from input_advanced
    npoints: int = 0
    n_init: int = 0
    dist_fact: float = 0
    flag_threshold: int = 0
    a_beta: float = 0
    b_beta: float = 0
    max_aspect_ratio: float = 0
    saveraster_flag: int = 0
    aspect_ratio_coeff: float = 0
    start_from_dist_flag: int = 0
    force_max_length: int = 0
    max_length: float = 0
    n_check_loop: int = 0
    restart_files: Optional[List[str]] = None
    restart_filling_parameters: Optional[List[float]] = None
    };

}