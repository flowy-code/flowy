#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include <filesystem>
#include <optional>

namespace Flowy
{

struct AscCrop
{
    double x_min;
    double x_max;
    double y_min;
    double y_max;
};

class AscFile
{
public:
    AscFile() = default;
    explicit AscFile( const std::filesystem::path & path, std::optional<AscCrop> crop = std::nullopt );
    void save( const std::filesystem::path & path );

    Vector2 lower_left_corner = { 0, 0 }; // Coordinates of lower left corner
    double cell_size          = 0;        // side length of square cell
    double no_data_value      = -9999;    // number that indicates lack of data

    // NOTE: that the order of rows in the asc file is opposite to the order of rows in the height_data
    // This is because we want the first row to correspond to the *low* y-values
    MatrixX height_data{}; // array that contains height data
    VectorX
        x_data{}; // one dimensional coordinates of the sampling grid in x direction (the lower left corner of each pixel)
    VectorX
        y_data{}; // one dimensional coordinates of the sampling grid in y direction (the lower left corner of each pixel)
};

} // namespace Flowy
