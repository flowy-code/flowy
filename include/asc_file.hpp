#pragma once
#include "definitions.hpp"

namespace Flowtastic
{
class AscFile
{
public:
    Vector2 lower_left_corner = { 0, 0 }; // Coordinates of lower left corner
    double cell_size          = 0;        // side length of square cell
    double no_data_value      = -9999;    // number that indicates lack of data
    MatrixX height_data{};                // array that contains height data
    VectorX x_data;                       // one dimensional coordinates of the sampling grid in x direction
    VectorX y_data;                       // one dimensional coordinates of the sampling grid in y direction
};
} // namespace Flowtastic