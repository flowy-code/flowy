#pragma once
#include "definitions.hpp"
#include "math.hpp"
#include <fmt/format.h>
#include <optional>
#include <vector>

namespace Flowtastic
{

class Lobe
{
private:
    double azimuthal_angle     = 0; // The azimuthal angle of the major semi-axis with respect to the x-axis
    double sin_azimuthal_angle = 0;
    double cos_azimuthal_angle = 1.0;

public:
    void set_azimuthal_angle( double azimuthal_angle )
    {
        this->azimuthal_angle = azimuthal_angle;
        sin_azimuthal_angle   = std::sin( azimuthal_angle );
        cos_azimuthal_angle   = std::cos( azimuthal_angle );
    }

    double get_azimuthal_angle() const
    {
        return azimuthal_angle;
    }

    double get_cos_azimuthal_angle() const
    {
        return cos_azimuthal_angle;
    }

    double get_sin_azimuthal_angle() const
    {
        return sin_azimuthal_angle;
    }

    Vector2 center                = { 0, 0 };     // The center of the ellipse
    Vector2 semi_axes             = { 1, 1 };     // The length of the semi-axies, first the major then the minor
    int dist_n_lobes              = 0;            // The distance to the initial lobe, counted in number of lobes
    int n_descendents             = 0;            // The cumulative number of descendent lobes
    std::optional<int> idx_parent = std::nullopt; // The idx of the parent lobe
    double alpha_inertial         = 0;            // The coefficient for the inertial contribution
    double thickness              = 0;            // Each lobe has a thickness

    // Checks if the point lies inside the lobe
    inline bool is_point_in_lobe( const Vector2 & coordinates ) const
    {
        // Translate the coordinates by the center of ellipse
        const double x = coordinates[0] - center[0];
        const double y = coordinates[1] - center[1];

        // Transform to the coordinate system centered on the center of the ellipse
        // To use eqn (x'/a)^2 + (y'/b)^2
        // The rotation matrix is [  {cos phi x, sin phi x}, {-sin phi x, cos phi x} ]
        const double x_prime = cos_azimuthal_angle * x + sin_azimuthal_angle * y;
        const double y_prime = -sin_azimuthal_angle * x + cos_azimuthal_angle * y;

        const double ellipse_lhs = ( x_prime / semi_axes[0] ) * ( x_prime / semi_axes[0] )
                                   + ( y_prime / semi_axes[1] ) * ( y_prime / semi_axes[1] );
        // Check if the point lies in the lobe
        return ellipse_lhs <= 1;
    }

    // Checks if a line segment between x1 and x2 intersects the lobe
    inline bool line_segment_intersects( const Vector2 & x1, const Vector2 & x2 ) const
    {
        // This matrix transforms coordinates into the axes frame of the ellipse
        // clang-format off
        const MatrixX ellipse_coordinate_transform = { { cos_azimuthal_angle, sin_azimuthal_angle },
                                                 { -sin_azimuthal_angle,cos_azimuthal_angle } };
        //clang-format on

        const Vector2 x1t = x1 - center;
        const Vector2 x2t = x2 - center;

        // Previously, this used xt::linalg::dot, but the handwritten version is much faster
        const Vector2 x1_prime = {ellipse_coordinate_transform(0,0) * x1t[0] + ellipse_coordinate_transform(0,1) * x1t[1] , ellipse_coordinate_transform(1,0) * x1t[0] + ellipse_coordinate_transform(1,1) * x1t[1]  };
        const Vector2 x2_prime = {ellipse_coordinate_transform(0,0) * x2t[0] + ellipse_coordinate_transform(0,1) * x2t[1] , ellipse_coordinate_transform(1,0) * x2t[0] + ellipse_coordinate_transform(1,1) * x2t[1]  };

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


    // Gives a point on the perimeter of the ellipse
    // The angle is relative to the semi major axis angle
    inline Vector2 point_at_angle(const double phi) const
    {
        const double a = semi_axes[0]; // major axis
        const double b = semi_axes[1]; // minor axis
        Vector2 coord = {a * std::cos(phi)*std::cos(azimuthal_angle) - b*std::sin(phi)*std::sin(azimuthal_angle), a*std::cos(phi)*std::sin(azimuthal_angle) + b*std::sin(phi)*std::cos(azimuthal_angle)};
        return coord + center;
    }

    inline std::vector<Vector2> rasterize_perimeter(int n_raster_points) const
    {
        auto phi_list = xt::linspace<double>(0.0, 2.0*Math::pi, n_raster_points, false);
        auto res = std::vector<Vector2>(n_raster_points);

        for (int idx_phi=0; idx_phi<n_raster_points; idx_phi++)
        {
            res[idx_phi] = point_at_angle(phi_list[idx_phi]);
        }

        return res;
    }
};

} // namespace Flowtastic