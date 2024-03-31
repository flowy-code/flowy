#include "topography.hpp"
#include "definitions.hpp"

namespace Flowtastic
{

Vector3 Topography::slope( const Vector2 & coordinates )
{
    // Have to find the four grid points that enclose the coordinates
    const int idx_x_lower  = int(coordinates[0] / cell_size);
    const int idx_x_higher = int(coordinates[0] / cell_size + 1);
    const int idx_y_lower  = int(coordinates[1] / cell_size);
    const int idx_y_higher = int(coordinates[1] / cell_size + 1);

    // We interpolate the height data with a function f(x,y) = alpha * x + beta * y + gamma * x * y + c
    // We want f(x,y) to coincide with the actual height data at the four enclosing corners
    // The height data is
    const double z00 = height_data(idx_x_lower, idx_y_lower);
    const double z01 = height_data(idx_x_higher, idx_y_lower);
    const double z10 = height_data(idx_x_lower, idx_y_higher);
    const double z11 = height_data(idx_x_higher, idx_y_higher);

    // And the values of x ad y are
    const double x0 = x_data[idx_x_lower];
    const double x1 = x_data[idx_x_higher];
    const double y0 = x_data[idx_y_lower];
    const double y1 = x_data[idx_y_higher];
    // We transform the coordinates into the intervala


    // This results in four equations, from which the values of alpha, beta and gamma can be derived
    // f(0,0) = z00    
    // f(1,0) = z10
    // f(0,1) = z01
    // f(1,1) = z11
    // alpha = f(1,0) - f(0,0)
    // beta = f(0,1) - f(0,0)
    // gamma = f(1,1) - f(1,0) + f(0,0) - f(0,1)


}

} // namespace Flowtastic