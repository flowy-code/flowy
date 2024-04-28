#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#ifdef WITH_NETCDF

#include "flowy/include/topography_file.hpp"
#include <optional>
#include <string>

namespace Flowy
{

class NetCDFFile : public TopographyFile
{
    std::optional<double> scale_factor = std::nullopt;
    std::optional<double> add_offset   = std::nullopt;
    void determine_scale_and_offset();
    MatrixX get_elevation_from_netcdf( int ncid );

public:
    bool compression      = false;
    bool shuffle          = true;
    int compression_level = 5;

    StorageDataType data_type = StorageDataType::Short;
    NetCDFFile()              = default;
    NetCDFFile( const Topography & topography, OutputQuantitiy output )
            : TopographyFile::TopographyFile( topography, output )
    {
    }

    explicit NetCDFFile(
        const std::filesystem::path & path, const std::optional<TopographyCrop> & crop = std::nullopt );

    std::string suffix() override
    {
        return ".nc";
    }

    void save( const std::filesystem::path & path ) override;
};

} // namespace Flowy

#endif
