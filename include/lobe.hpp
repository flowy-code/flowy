#pragma once
#include "definitions.hpp"
#include <fmt/format.h>
#include <optional>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#include <xtensor-blas/xlinalg.hpp>
#pragma GCC diagnostic pop

namespace Flowtastic
{

class Lobe
{
public:
    double azimuthal_angle        = 0;        // The azimuthal angle of the major semi-axis with respect to the x-axis
    Vector2 center                = { 0, 0 }; // The center of the ellipse
    Vector2 semi_axes             = { 1, 1 }; // The length of the semi-axies, first the major then the minor
    int dist_n_lobes              = 0;        // The distance to the initial lobe, counted in number of lobes
    int n_descendents             = 0;        // The cumulative number of descendent lobes
    std::optional<int> idx_parent = std::nullopt; // The idx of the parent lobe
    double alpha_inertial         = 0;            // The coefficient for the inertial contribution
    double thickness              = 0;            // Each lobe has a thickness

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

    // Checks if a line segment between x1 and x2 intersects the lobe
    inline bool line_segment_intersects( const Vector2 & x1, const Vector2 & x2 ) const
    {
        // This matrix transforms coordinates into the axes frame of the ellipse
        // clang-format off
        MatrixX ellipse_coordinate_transform = { { std::cos( azimuthal_angle ), std::sin( azimuthal_angle ) },
                                                 { -std::sin( azimuthal_angle ),std::cos( azimuthal_angle ) } };
        //clang-format on

        const Vector2 x1_prime = xt::linalg::dot( ellipse_coordinate_transform, x1 - center  );
        const Vector2 x2_prime = xt::linalg::dot( ellipse_coordinate_transform, x2 - center  );
        const Vector2 diff     = x2_prime - x1_prime;

        const double a = semi_axes[0];
        const double b = semi_axes[1];
        // In the ellipse coordinates the equation for the perimeter is
        // (x/a)^2 + (y/b)^2 = 1

        // The line segment can be described by the vector equation
        // v(t) = x1 + diff*t for t in [0,1]

        // Plugging this into the ellipse equation, yields an equation for t
        // alpha * t^2 + beta * t + gamma = 0
        const double alpha = 1.0 / ( a * a ) * ( diff[0] * diff[0] ) + 1.0 / ( b * b ) * ( diff[1] * diff[1] );
        const double beta  = 2.0 * ( x1_prime[0] * diff[0] / ( a * a ) + x1_prime[1] * diff[1] / ( b * b ) );
        const double gamma = x1_prime[0] * x1_prime[0] / ( a * a ) + x1_prime[1] * x1_prime[1] / ( b * b ) - 1;

        // The solution to this quadratic equation is
        // t = (-beta +- sqrt(beta^2 - 4*alpha*gamma)) / (2*alpha)
        
        // Therefore, if beta^2 - 4*alpha*gamma < 0, the line segment misses the ellipse
        const double radicand = beta*beta - 4*alpha*gamma;
        if (radicand < 0)
            return false;

        // Else, we compute the two points of intersection and check if they fall into the interval [0, 1]
        const double t1 = (-beta - std::sqrt(radicand)) / (2.0*alpha);
        const double t2 = (-beta + std::sqrt(radicand)) / (2.0*alpha);

        if (t1 >= 0 && t1 <= 1.0)
            return true;

        if (t2 >= 0 && t2 <= 1.0)
            return true;

        return false;
    }
};

} // namespace Flowtastic