#pragma once
#include "definitions.hpp"

namespace Flowtastic
{

class Lobe
{
public:
    double azimuthal_angle = 0;        // The azimuthal angle of the major semi-axis with respect to the x-axis
    Vector2 center         = { 0, 0 }; // The center of the ellipse
    Vector2 semi_axes      = { 1, 1 }; // The length of the semi-axies, first the major then the minor
    int dist_n_lobes       = 0;        // The distance to the initial lobe, counted in number of lobes
    int n_descendents      = 0;        // The cumulative number of descendent lobes
    int idx_parent         = 0;        // The idx of the parent lobe
    double alpha_inertial  = 0;        // The coefficient for the inertial contribution
};

} // namespace Flowtastic