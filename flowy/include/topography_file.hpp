#pragma once
#include "flowy/include/definitions.hpp"
#include "flowy/include/topography.hpp"
#include <filesystem>
namespace Flowy
{

struct TopographyCrop
{
    double x_min;
    double x_max;
    double y_min;
    double y_max;
};

class TopographyFile
{
public:
    // NOTE: that the order of rows in the asc file is opposite to the order of rows in the height_data
    // This is because we want the first row to correspond to the *low* y-values
    MatrixX height_data{}; // array that contains height data
    VectorX
        x_data{}; // one dimensional coordinates of the sampling grid in x direction (the lower left corner of each pixel)
    VectorX
        y_data{}; // one dimensional coordinates of the sampling grid in y direction (the lower left corner of each pixel)
    double no_data_value = -9999; // number that indicates lack of data

    // Constructor that takes topography
    TopographyFile( const Topography & topography )
            : height_data( topography.height_data ), x_data( topography.x_data ), y_data( topography.y_data )
    {
    }

    // Cell size of the rectilinear grid
    double cell_size()
    {
        if( x_data.size() < 2 )
        {
            return 1.0;
        }
        return x_data[1] - x_data[0];
    }

    // Coordinates of lower left corner
    Vector2 lower_left_corner()
    {
        return { x_data[0], y_data[0] };
    }

    virtual Topography to_topography()                      = 0;
    virtual void save( const std::filesystem::path & path ) = 0;
    virtual ~TopographyFile()                               = default;

protected:
    // Crops the topography according to user-defined crop dimensions
    virtual void crop_topography( const TopographyCrop & crop_dims )
    {
        Vector2 old_lower_left_corner = lower_left_corner();
        // If cropping is used, we slice the height_data array
        int idx_x_min = std::clamp<int>(
            ( crop_dims.x_min - old_lower_left_corner[0] ) / cell_size(), 0, height_data.shape()[0] - 1 );
        int idx_x_max = std::clamp<int>(
            ( crop_dims.x_max - old_lower_left_corner[0] ) / cell_size(), 0, height_data.shape()[0] - 1 );
        int idx_y_min = std::clamp<int>(
            ( crop_dims.y_min - old_lower_left_corner[1] ) / cell_size(), 0, height_data.shape()[1] - 1 );
        int idx_y_max = std::clamp<int>(
            ( crop_dims.y_max - old_lower_left_corner[1] ) / cell_size(), 0, height_data.shape()[1] - 1 );

        height_data
            = xt::view( height_data, xt::range( idx_x_min, idx_x_max + 1 ), xt::range( idx_y_min, idx_y_max + 1 ) );

        // You need to modify x_data and y_data, since they are cropped now.
        x_data = xt::view( x_data, xt::range( idx_x_min, idx_x_max + 1 ) );
        y_data = xt::view( y_data, xt::range( idx_y_min, idx_y_max + 1 ) );
    }
};

} // namespace Flowy