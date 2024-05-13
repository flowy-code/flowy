#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include "flowy/include/netcdf_file.hpp"
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace Flowy::Config
{

struct OutputSettings
{
    bool crop_to_content      = false;
    bool use_netcdf           = false;
    bool compression          = true;
    int compression_level     = 5;
    bool shuffle              = true;
    StorageDataType data_type = StorageDataType::Float;
};

class InputParams
{
public:
    // ===================================================================================
    // new settings
    // ===================================================================================

    // The folder output is written to
    std::filesystem::path output_folder = "./output";

    OutputSettings output_settings{};

    // If set to true one csv file, per flow, is written to the output folder.
    // The files are named 'lobes_{idx_flow}.csv' and contain information about the lobes in that specific flow
    bool write_lobes_csv      = false;
    bool print_remaining_time = false;
    bool save_final_dem       = false;

    // The tolerance in the volume ratio when finding the threshold thickness for masking
    double masking_tolerance = 1e-5;
    // The maximum number of bisection search iterations when finding the threshold thickness for masking
    int masking_max_iter = 20;

    std::optional<int> rng_seed = std::nullopt;

    // Whether to apply a volume correction or not. By default this is false
    bool volume_correction = false;

    // ===================================================================================
    // MrLavaLoba settings, for feature parity
    // ===================================================================================

    std::string run_name{};                  // Name of the run (used to save the parameters and the output)
    std::filesystem::path source{};          // File name of ASCII digital elevation model (.asc file)
    std::vector<Vector2> vent_coordinates{}; // of shape [n_vents, 2]

    int n_vents() const
    {
        return vent_coordinates.size();
    }

    bool save_hazard_data{}; // If true, a raster map is saved, such that the values represent the probability of a cell
                             // to be covered by lava.
    int n_flows{ 1 };        // Number of flows per run
    int n_lobes{};           // Number of lobes per flow
    double thickening_parameter{}; // Parameter that affects calculation of the slope [0,1]
    std::optional<double>
        prescribed_lobe_area{}; // User defined (constant) lobe area. Either the lobe area or the average lobe thickness
                                // is computed, depending on the fixed_dimension_flag.
    std::optional<double> prescribed_avg_lobe_thickness{}; // User defined average lobe thickness

    // Variables we don't understand
    std::vector<double> masking_threshold{};
    int min_n_lobes{ 0 };
    int max_n_lobes{ 1 };

    /*
    Inertial exponent, used for calculating the inertial modification to the azimuthal angle (also depends on the slope,
    which cannot be negatives):
    inertial_exponent = 0: no inertial contribution is added
    direction; inertial_exponent > 0: some effect of the parent lobe semi-major axis direction is added to the budding
    lobe semi-major axis direction.
    */
    double inertial_exponent{};

    /*
    The lobe_exponent dictates the selection of the parent lobe. The closer that the `lobe_exponent` parameter is to 0, the
    larger is the probability that the new lobe will be generated from a
    newer lobe (which is some "distance", in number of lobes, away from the vent).
    lobe_exponent = 1: the parent lobe is chosen
    with a uniform probability distribution.
    lobe_exponent = 0: the new lobe is always generated from the last one.
    lobe_exponent between 0 and 1: a factor is chosen from a probability distribution given by P(x) = (1/eps * x^{1/eps
    -1 }), where eps is the lobe exponent
    */
    double lobe_exponent{};

    /*
    `max_slope_prob` dictates the extent of the azimuthal angle perturbation.
    When max_slope_prob = 0: the angle perturbation is drawn from a uniform probability distribution
    max_slope_prob > 0: the angle perturbation is drawn from a truncated normal distribution, constrained within
    [-pi,pi] max_slope_prob = 1: the angle is not perturbed.
    */
    double max_slope_prob{};

    /*
    thickness_ratio =  (thickness of the first lobe)/ (thickness of the last lobe)
    thickness_ratio < 1   => the thickness increases with lobe "age"
    thickness_ratio = 1   => all the lobes have the same thickness
    thickness_ratio > 1   => the thickness decreases with lobe "age"
    */
    double thickness_ratio{};

    /*
    This flag selects which dimension of the lobe is fixed:
    fixed_dimension_flag = 1  => the area of the lobes is assigned
    fixed_dimension_flag = 2  => the thickness of the lobes is assigend
    */
    int fixed_dimension_flag{};

    /*
    See VentFlag class
    */
    int vent_flag{};

    /*To define the ends of fissures */
    std::optional<std::vector<Vector2>> fissure_end_coordinates{}; // of shape [n_vents, 2]
    /*Fissure probabilities, should be of length n_vents (or number of fissures)*/
    std::optional<std::vector<double>> fissure_probabilities{};
    std::optional<double> total_volume{};
    std::optional<double> east_to_vent{};
    std::optional<double> west_to_vent{};
    std::optional<double> south_to_vent{};
    std::optional<double> north_to_vent{};
    std::optional<std::string> channel_file{};
    std::optional<double> alfa_channel{};
    std::optional<double> d1{};
    std::optional<double> d2{};
    std::optional<double> eps{};
    std::optional<std::string> union_diff_file{};

    // ===================================================================================
    // Advanced MrLavaLoba settings, for feature parity
    // ===================================================================================

    int npoints{ 30 }; // Number of points for rasterizing the ellipse
    int n_init{ 0 };   // Number of repetitions of the first lobe (useful for initial spreading)

    /*This factor is to choose where the center of the new lobe will be:
      dist_fact = 0 : the center of the new budding lobe is on the budding point, which is on the perimeter of the
      parent lobe; dist fact > 0 : this then increases the distance of the center of the new lobe from the budding point
      dist_fact = 1 : the parent lobe and budding lobe touch in one point only (the budding points)
      dist_fact > 1 : the parent lobe and budding lobe do not touch
      By default the parent lobe and budding lobe overlap such that they intersect at a point (the budding point) */
    double dist_fact{ 1 };

    /*
    Not used/implemented?
    */
    int flag_threshold{};

    /*
    We haven't implemented this, but the number of lobes could be drawn from a beta probability density function (not a
    beta probability distribution)
    */
    double a_beta{ 0 };
    double b_beta{ 0 };

    double max_aspect_ratio{}; // Maximum possible aspect ratio of the lobes
    int saveraster_flag{};     // We don't use this

    /*
    Decides how the lobe aspect ratio scales with the slope
    aspect_ratio_coeff = 0: the lobe is always a circle because the slope does not affect the aspect ratio at all
    aspect_ratio_coeff > 0: the lobe is an ellipse, such that the aspect ratio increases with the slope; therefore the
    lobes become more elongated with steeper slopes
    */
    double aspect_ratio_coeff{};

    /*
    Not implemented
    */
    int start_from_dist_flag{};

    int force_max_length{}; // Not implemented

    /*
    Not implemented
    */
    double max_length{};

    /*
    Not implemented
    */
    int n_check_loop{};

    /*
    This is a list of names of existing *_thickness_*.asc files, to be used if
    you want to start a simulation from previous
    flows.
    */
    std::optional<std::vector<std::filesystem::path>> restart_files{};

    // Multiply restart file thickness by a constant value
    std::optional<std::vector<double>> restart_filling_parameters{};
};

} // namespace Flowy::Config
