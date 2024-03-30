#pragma once
#include "definitions.hpp"
#include <filesystem>
#include <optional>
#include <string>

namespace Flowtastic::Config
{

class InputParams
{
public:
    std::filesystem::path run_name{}; // Name of the run (used to save the parameters and the output)
    std::filesystem::path source{};   // File name of ASCII digital elevation model (.asc file)
    MatrixX vent_coordinates{};       // of shape [n_vents, 2]
    bool save_hazard_data{}; // If true, a raster map is saved, such that the values represent the probability of a cell
                             // to be covered by lava.
    int n_flows{};           // Number of flows per run
    int n_lobes{};           // Number of lobes per flow
    double thickening_parameter{}; // Parameter that affects calculation of the slope [0,1]
    std::optional<double>
        prescribed_lobe_area{}; // User defined (constant) lobe area. Either the lobe area or the average lobe thickness
                                // is computed, depending on the fixed_dimension_flag.
    std::optional<double> prescribed_avg_lobe_thickness{}; // User defined average lobe thickness

    // Variables we don't understand
    // masking_threshold: float = 0
    int min_n_lobes{};
    int max_n_lobes{};

    /*
    Inertial exponent, used for calculating the inertial modification to the azimuthal angle (also depends on the slope,
    which cannot be negatives): inertial_exponent = 0 => the max probability direction for the new lobe is the max slope
    direction; inertial_exponent > 0 => the max probability direction for the new lobe takes into account also the
    direction of the parent lobe and the inertia increaes with increasing exponent
    */
    float inertial_exponent{};

    /*
    The lobe_exponent is associated with the probability that a new lobe will
    be descended from a "young" or "old" parent lobe (far away from or close to the vent, when the
    flag `start_from_dist_flag` is set). The closer that `lobe_exponent` is to 0, the
    larger is the probability that the new lobe will be generated from a
    "younger" lobe (which is some "distance", in number of lobes,s away from the vent).
    lobe_exponent = 1 => th                                                                  e parent lobe is chosen
    with a uniform probability distribution. lobe_exponent = 0 => the new lobe is generated from the last one.
    */
    float lobe_exponent{};

    /*
    `max_slope_prob` is related to the probability that the direction of
    the new lobe is close to the maximum slope direction:
    max_slope_prob = 0 => all the directions have the same probability;
    max_slope_prob > 0 => the maximum slope direction has a larger
                          probability, and the probability increases with increasing
                          value of the parameter;
    max_slope_prob = 1 => the direction of the new lobe is the maximum slope direction.
    */
    float max_slope_prob{};

    /*
    thickness_ratio =  (thickness of the first lobe)/ (thickness of the last lobe)
    thickness_ratio < 1   => the thickness increases with lobe "age"
    thickness_ratio = 1   => all the lobes have the same thickness
    thickness_ratio > 1   => the thickness decreases with lobe "age"
    */
    float thickness_ratio{};

    /*
    This flag selects which dimension of the lobe is fixed:
    fixed_dimension_flag = 1  => the area of the lobes is assigned
    fixed_dimension_flag = 2  => the thickness of the lobes is assigend
    */
    int fixed_dimension_flag{};

    /*
    This flag select how multiple initial coordinates are treated:
    vent_flag = 0 => the initial lobes are on the vents coordinates
                      and the flows start initially from the first vent,
                      then from the second and so on.
    vent_flag = 1 => the initial lobes are chosen randomly from the vents
                      coordinates and each vent has the same probability
    vent_flag = 2 => the initial lobes are on the polyline connecting
                      the vents and all the point of the polyline
                      have the same probability
    vent_flag = 3 => the initial lobes are on the polyline connecting
                      the vents and all the segments of the polyline
                      have the same probability
    vent_flag = 4 => the initial lobes are on multiple
                      fissures and all the point of the fissures
                      have the same probability
    vent_flag = 5 => the initial lobes are on multiple
                      fissures and all the fissures
                      have the same probability
    vent_flag = 6 => the initial lobes are on the polyline connecting
                      the vents and the probability of
                      each segment is fixed by "fissure probabilities"
    vent_flag = 7 => the initial lobes are on multiple
                      fissures and the probability of
                      each fissure is fixed by "fissure_probabilities"
    vent_flag = 8 => the initial lobes are chosen randomly from the vents
                      coordinates and the probability of each vent
    */
    int vent_flag{};

    std::optional<float> fissure_probabilities{};
    std::optional<float> total_volume{};
    std::optional<float> east_to_vent{};
    std::optional<float> west_to_vent{};
    std::optional<float> south_to_vent{};
    std::optional<float> north_to_vent{};
    std::optional<std::string> channel_file{};
    std::optional<float> alfa_channel{};
    std::optional<float> d1{};
    std::optional<float> d2{};
    std::optional<float> eps{};
    std::optional<std::string> union_diff_file{};

    // from input_advanced
    int npoints{}; // Number of points for the lobe profile
    int n_init{};  // Number of repetitions of the first lobe (useful for initial spreading)

    /*This factor is to choose where the center of the new lobe will be:
      dist_fact = 0 => the center of the new lobe is on the border of the
                        previous one;
      dist fact > 0 => increase the distance of the center of the new lobe
                        from the border of the previous one;
      dist_fact = 1 => the two lobes touch in one point only.*/
    float dist_fact{};

    /*
    Flag to select if it is cutted the volume of the area
    flag_threshold = 1  => volume
    flag_threshold = 2  => area
    */
    int flag_threshold{};

    /*
    The number of lobes of the flow is defined accordingly to a random uniform
    distribution or to a beta law, as a function of the flow number.
    a_beta, b_beta = 0 => n_lobes is sampled randomly in [min_n_lobes,max_n_lobes]
    a_beta, b_beta > 0 => n_lobes = min_n_lobes + 0.5 * ( max_n_lobes - min_n_lobes )
                                                 * beta(flow/n_flows,a_beta,b_beta)
    */
    float a_beta{};
    float b_beta{};

    float max_aspect_ratio{}; // Maximum aspect ration of the lobes
    int saveraster_flag{};    // if saveraster_flag = 1 then the raster output is saved as a *.asc file

    /*
    This parameter affect the shape of the lobes. The larger is this parameter
    the larger is the effect of a small slope on the eccentricity of the lobes:
    aspect_ratio_coeff = 0 => the lobe is always a circle
    aspect_ratio_coeff > 0 => the lobe is an ellipse, with the aspect ratio
                              increasing with the slope
    */
    float aspect_ratio_coeff{};

    /*
    This flag controls which lobes have larger probability:
    start_from_dist_flag = 1 => the lobes with a larger distance from
                                 the vent have a higher probability
    start_form_dist_flag = 0 => the younger lobes have a higher
                                 probability
    */
    int start_from_dist_flag{};

    int force_max_length{}; // Flag for maximum distances (number of chained lobes) from the vent

    /*
    Maximum distances (number of chained lobes) from the vent
    This parameter is used only when force_max_length = 1
    */
    float max_length{};

    /*
    This parameter is to avoid a chain of loop getting stuck in a hole. It is
    active only when n_check_loop>0. When it is greater than zero the code
    check if the last n_check_loop lobes are all in a small box. If this is the
    case then the slope modified by flow is evaluated, and the hole is filled.
    */
    int n_check_loop{};

    /*
    This is a list of names of existing *_thickness_*.asc files, to be used if
    you want to start a simulation from previous
    flows.
    */
    std::optional<std::vector<std::filesystem::path>> restart_files{};

    // This seems to be described nowhere
    std::optional<std::vector<std::filesystem::path>> restart_filling_parameters{};
};

} // namespace Flowtastic::Config