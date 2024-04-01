#pragma once
#include "definitions.hpp"
#include "fmt/ostream.h"
#include <fmt/format.h>

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

    // Checks if the point lies inside the lobe
    inline bool is_point_in_lobe( const Vector2 & coordinates ) const
    {
        // Translate the coordinates by the center of ellipse
        double x = coordinates[0] - center[0];
        double y = coordinates[1] - center[1];

        // Transform to the coordinate system centered on the center of the ellipse
        // To use eqn (x'/a)^2 + (y'/b)^2
        // The rotation matrix is [  {cos phi x, sin phi x}, {-sin phi x, cos phi x} ]
        double x_prime = std::cos( azimuthal_angle ) * x + std::sin( azimuthal_angle ) * y;
        double y_prime = -std::sin( azimuthal_angle ) * x + std::cos( azimuthal_angle ) * y;

        double ellipse_lhs = ( x_prime / semi_axes[0] ) * ( x_prime / semi_axes[0] )
                             + ( y_prime / semi_axes[1] ) * ( y_prime / semi_axes[1] );
        // Check if the point lies in the lobe
        return ellipse_lhs <= 1;
    }

    // Determines the rectangular "bounding box" enclosing the ellipse (which can be tilted).
    // Returns the coordinates of the vertices of the bounding box
    // clockwise around the rectangle
    inline std::vector<Vector2> bounding_box() const
    {
        Vector2 major_axis_vec
            = { semi_axes[0] * std::cos( azimuthal_angle ), semi_axes[0] * std::sin( azimuthal_angle ) };
        Vector2 minor_axis_vec
            = { -semi_axes[1] * std::sin( azimuthal_angle ), semi_axes[1] * std::cos( azimuthal_angle ) };

        // The points are given by center plus/minus major_axis_vec plus/minus minor_axis_vec
        Vector2 point1
            = { center[0] + major_axis_vec[0] - minor_axis_vec[0], center[1] + major_axis_vec[1] - minor_axis_vec[1] };
        Vector2 point2
            = { center[0] + major_axis_vec[0] + minor_axis_vec[0], center[1] + major_axis_vec[1] + minor_axis_vec[1] };
        Vector2 point3
            = { center[0] - major_axis_vec[0] + minor_axis_vec[0], center[1] - major_axis_vec[1] + minor_axis_vec[1] };
        Vector2 point4
            = { center[0] - major_axis_vec[0] - minor_axis_vec[0], center[1] - major_axis_vec[1] - minor_axis_vec[1] };

        return { point1, point2, point3, point4 };
    }
};

} // namespace Flowtastic