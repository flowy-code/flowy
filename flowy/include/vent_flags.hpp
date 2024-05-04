#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include "flowy/include/math.hpp"
#include <fmt/format.h>
#include <optional>
#include <random>
#include <stdexcept>
#include <vector>

namespace Flowy
{
/*
Abstract class for handling different input.vent_flag options.
By default, the fissure chosen is selected from a uniform probability distribution.
Therefore, all fissures have the same probability.
The lobe position is chosen using the fissure end points (using a uniform probability distribution),
but get_fissure_endpoints needs to be defined in the child VentFlag class if fissures need to be calculated.
*/
class VentFlag
{
public:
    using fissProbT  = std::vector<double>;
    using ventCoordT = std::vector<Vector2>;

    VentFlag(
        const ventCoordT & vent_coordinates, const std::optional<fissProbT> & fissure_probabilities,
        const std::optional<ventCoordT> & fissure_end_coordinates, std::mt19937 & gen )
            : vent_coordinates( vent_coordinates ),
              fissure_probabilities( fissure_probabilities ),
              fissure_end_coordinates( fissure_end_coordinates ),
              gen( gen )
    {
    }

    int n_vents()
    {
        return vent_coordinates.size();
    }

    virtual int n_fissures()
    {
        return 0;
    }

    virtual void compute_line_segment_weights()
    {
        fissure_weights = std::vector<double>( n_fissures(), 1.0 );
    }

    virtual void get_fissure_endpoints( int idx_fissure [[maybe_unused]] ) {}

    // Get the lobe position from the endpoints of the fissure
    [[nodiscard]] virtual Vector2 get_position_from_endpoints()
    {
        return x1 + ( 1.0 - dist_line( gen ) ) * ( x2 - x1 );
    }

    virtual void validate_fissure_probabilities()
    {
        if( fissure_probabilities.has_value() )
        {
            if( fissure_probabilities->size() != n_fissures() )
            {
                throw std::runtime_error( fmt::format(
                    "The size of fissure_probabilities (= {}), does not match the number of fissures (= {})",
                    fissure_probabilities->size(), n_fissures() ) );
            }
        }
        else
        {
            throw std::runtime_error( "Fissure probabilities need to be set!" );
        }
    }

    virtual void validate_fissure_end_coords()
    {
        if( !fissure_end_coordinates.has_value() )
        {
            throw std::runtime_error( "x_vent_end and y_vent_end need to be set!" );
        }
    }

    [[nodiscard]] virtual Vector2 get_position()
    {
        compute_line_segment_weights();
        std::discrete_distribution<int> dist_segment( fissure_weights.begin(), fissure_weights.end() );
        int idx_fissure = dist_segment( gen );
        get_fissure_endpoints( idx_fissure );
        return get_position_from_endpoints();
    }

    virtual ~VentFlag() = default;

protected:
    ventCoordT vent_coordinates{};
    std::optional<fissProbT> fissure_probabilities{};
    std::optional<ventCoordT> fissure_end_coordinates{};
    std::vector<double> fissure_weights{};

    Vector2 x1{};
    Vector2 x2{};

    std::uniform_real_distribution<double> dist_line = std::uniform_real_distribution<double>( 0, 1 );
    std::mt19937 & gen;
};

/*
The initial lobes are on the vents coordinates and the flows start initially from the first vent,
then from the second and so on.
*/
class VentFlag0 : public VentFlag
{
    int idx_flow{};
    int n_flows{};

public:
    VentFlag0( int idx_flow, int n_flows, const ventCoordT & vent_coordinates, std::mt19937 & gen )
            : VentFlag( vent_coordinates, std::nullopt, std::nullopt, gen ), idx_flow( idx_flow ), n_flows( n_flows )
    {
    }

    Vector2 get_position() override
    {
        int idx_vent = std::floor( idx_flow * n_vents() / n_flows );
        return vent_coordinates[idx_vent];
    }
};

/*
The initial lobes are chosen randomly from the vent coordinates, and each vent has the same probability
*/
class VentFlag1 : public VentFlag
{
public:
    using VentFlag::VentFlag;

    Vector2 get_position() override
    {
        auto dist    = std::uniform_int_distribution<int>( 0, n_vents() - 1 );
        int idx_vent = dist( gen );
        return vent_coordinates[idx_vent];
    }
};

/*
The initial lobes are on the polyline connecting the vents.
All the segments of the polyline have the same probability (in other words, each "fissure" is a line seg)
*/
class VentFlag3 : public VentFlag
{
public:
    using VentFlag::VentFlag;

    int n_fissures() override
    {
        return n_vents() - 1;
    }

    void get_fissure_endpoints( int idx_fissure ) override
    {
        x1 = vent_coordinates[idx_fissure];
        x2 = vent_coordinates[idx_fissure + 1];
    }
};

/*
The initial lobes are on the polyline connecting the vents.
All fissures (line segments connecting the vents) have a probability that is weighted by the length of the "fissure".
Once a fissure has been chosen, all points in the fissure have the same probability, and the location of the lobe is
selected using a uniform probability distribution
*/
class VentFlag2 : public VentFlag3
{
public:
    using VentFlag3::VentFlag3;

    void compute_line_segment_weights() override
    {
        fissure_weights = std::vector<double>( n_fissures() );
        for( int i_fissure = 0; i_fissure < n_fissures(); i_fissure++ )
        {
            fissure_weights[i_fissure]
                = xt::linalg::norm( vent_coordinates[i_fissure + 1] - vent_coordinates[i_fissure] );
        }
    }
};

/*
The initial lobes are on multiple fissures and all the fissures have the same probability.
Therefore the fissure is chosen by using a unifrom probability distribution.
The lobe position is chosen using the fissure end points (using a uniform probability distribution)
*/
class VentFlag5 : public VentFlag
{
public:
    using VentFlag::VentFlag;

    int n_fissures() override
    {
        return n_vents();
    }

    void get_fissure_endpoints( int idx_fissure ) override
    {
        validate_fissure_end_coords();
        x1 = vent_coordinates[idx_fissure];
        x2 = fissure_end_coordinates.value()[idx_fissure];
    }
};

/*
The initial lobes are on multiple fissures, but the probability of choosing a fissure is weighted by its length.
The lobe position is chosen using the fissure end points (using a uniform probability distribution)
*/
class VentFlag4 : public VentFlag5
{
public:
    using VentFlag5::VentFlag5;

    void compute_line_segment_weights() override
    {
        validate_fissure_end_coords();

        fissure_weights = std::vector<double>( n_fissures() );
        for( int i_fissure = 0; i_fissure < n_fissures(); i_fissure++ )
        {
            fissure_weights[i_fissure]
                = xt::linalg::norm( fissure_end_coordinates.value()[i_fissure] - vent_coordinates[i_fissure] );
        }
    }
};

/*
The initial lobes are on the polyline connecting the vents (that's our "fissure")
and the probability of each segment is fixed by "fissure probabilities"
*/
class VentFlag6 : public VentFlag3
{
public:
    using VentFlag3::VentFlag3;

    void compute_line_segment_weights() override
    {
        validate_fissure_probabilities();
        fissure_weights = fissure_probabilities.value();
    }
};

/*
The initial lobes are on multiple fissures and the probability of each fissure is fixed by "fissure_probabilities"
The lobe position is chosen using the fissure end points (using a uniform probability distribution)
*/
class VentFlag7 : public VentFlag5
{
public:
    using VentFlag5::VentFlag5;

    void compute_line_segment_weights() override
    {
        validate_fissure_probabilities();
        fissure_weights = fissure_probabilities.value();
    }
};

/*
The initial lobes are chosen randomly from the vent coordinates, as well as from the fissure probabilities
*/
class VentFlag8 : public VentFlag
{
public:
    using VentFlag::VentFlag;
    int n_fissures() override
    {
        return n_vents();
    }

    Vector2 get_position() override
    {
        validate_fissure_probabilities();
        std::discrete_distribution<int> dist_vent( fissure_probabilities->begin(), fissure_probabilities->end() );
        int idx_vent = dist_vent( gen );
        return vent_coordinates[idx_vent];
    }
};

} // namespace Flowy
