#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include "flowy/include/topography.hpp"
#include <fmt/ostream.h>
#include <cstddef>
#include <filesystem>
#include <stdexcept>
#include <string>
namespace Flowy
{

struct TopographyCrop
{
    double x_min;
    double x_max;
    double y_min;
    double y_max;
};

enum class OutputQuantitiy
{
    Hazard,
    Height
};

class TopographyFile
{
protected:
    virtual std::string suffix() = 0;

    std::filesystem::path handle_suffix( const std::filesystem::path & path_in )
    {
        std::filesystem::path path_copy = path_in;
        auto filename_old               = path_in.filename();
        return path_copy.replace_filename( fmt::format( "{}{}", filename_old.string(), suffix() ) );
    }

public:
    // NOTE: that the order of rows in the asc file is opposite to the order of rows in the height_data
    // This is because we want the first row to correspond to the *low* y-values
    VectorX
        x_data{}; // one dimensional coordinates of the sampling grid in x direction (the lower left corner of each pixel)
    VectorX
        y_data{}; // one dimensional coordinates of the sampling grid in y direction (the lower left corner of each pixel)
    MatrixX data{}; // array that contains height data

    double no_data_value = -9999; // number that indicates lack of data

    TopographyFile() = default;

    // Constructor that takes topography
    TopographyFile( const Topography & topography, OutputQuantitiy output )
            : x_data( topography.x_data ), y_data( topography.y_data ), no_data_value( topography.no_data_value )
    {
        if( output == OutputQuantitiy::Height )
        {
            data = topography.height_data;
        }
        else if( output == OutputQuantitiy::Hazard )
        {
            data = topography.hazard;
        }
        else
        {
            throw std::runtime_error( "Unknown data type" );
        }
    }

    // Cell size of the rectilinear grid
    double cell_size() const
    {
        if( x_data.size() < 2 )
        {
            return 1.0;
        }
        return std::abs( x_data[1] - x_data[0] );
    }

    // Coordinates of lower left corner
    Vector2 lower_left_corner() const
    {
        if( x_data.size() < 1 || y_data.size() < 1 )
        {
            return { 0, 0 };
        }

        return { x_data[0], y_data[0] };
    }

    Topography to_topography() const
    {
        return Topography( data, x_data, y_data, no_data_value );
    }

    virtual void save( const std::filesystem::path & path ) = 0;

    virtual ~TopographyFile() = default;

    // Crops the topography according to user-defined crop dimensions
    virtual void crop_to_content()
    {
        size_t idx_x_min = data.shape()[0] - 1;
        size_t idx_y_min = data.shape()[1] - 1;
        size_t idx_x_max{};
        size_t idx_y_max{};

        for( size_t ix = 0; ix < data.shape()[0]; ix++ )
        {
            for( size_t iy = 0; iy < data.shape()[1]; iy++ )
            {
                auto d = data( ix, iy );

                if( d > no_data_value + 1e-3 )
                {
                    if( ix < idx_x_min )
                        idx_x_min = ix;
                    if( ix > idx_x_max )
                        idx_x_max = ix;
                    if( iy < idx_y_min )
                        idx_y_min = iy;
                    if( iy > idx_y_max )
                        idx_y_max = iy;
                }
            }
        }

        // TODO: figure out what to do if idx_min == idx_max
        if( idx_x_min == idx_x_max || idx_y_min == idx_y_max )
        {
            return;
        }

        data   = xt::view( data, xt::range( idx_x_min, idx_x_max + 1 ), xt::range( idx_y_min, idx_y_max + 1 ) );
        x_data = xt::view( x_data, xt::range( idx_x_min, idx_x_max + 1 ) );
        y_data = xt::view( y_data, xt::range( idx_y_min, idx_y_max + 1 ) );
    }

    std::string to_string()
    {
        return fmt::format(
            "no_data_value = {}\nx_data = {}\ny_data = {}\nll_corner = {}\ncell_size = {}\ndata = {}\n", no_data_value,
            fmt::streamed( x_data ), fmt::streamed( y_data ), fmt::streamed( lower_left_corner() ), cell_size(),
            fmt::streamed( data ) );
    }

protected:
    // Crops the topography according to user-defined crop dimensions
    virtual void crop_topography( const TopographyCrop & crop_dims )
    {
        Vector2 old_lower_left_corner = lower_left_corner();
        // If cropping is used, we slice the data array
        int idx_x_min
            = std::clamp<int>( ( crop_dims.x_min - old_lower_left_corner[0] ) / cell_size(), 0, data.shape()[0] - 1 );
        int idx_x_max = std::clamp<int>(
            ( crop_dims.x_max - old_lower_left_corner[0] ) / cell_size() + 1, 0, data.shape()[0] - 1 );
        int idx_y_min
            = std::clamp<int>( ( crop_dims.y_min - old_lower_left_corner[1] ) / cell_size(), 0, data.shape()[1] - 1 );
        int idx_y_max = std::clamp<int>(
            ( crop_dims.y_max - old_lower_left_corner[1] ) / cell_size() + 1, 0, data.shape()[1] - 1 );

        data = xt::view( data, xt::range( idx_x_min, idx_x_max + 1 ), xt::range( idx_y_min, idx_y_max + 1 ) );

        // You need to modify x_data and y_data, since they are cropped now.
        x_data = xt::view( x_data, xt::range( idx_x_min, idx_x_max + 1 ) );
        y_data = xt::view( y_data, xt::range( idx_y_min, idx_y_max + 1 ) );
    }
};

} // namespace Flowy
