#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers

#include "flowy/include/topography_file.hpp"
#include <optional>

namespace Flowy
{

enum StorageDataType
{
    Short,
    Float,
    Double
};

class NetCDFFile : public TopographyFile
{
    std::optional<double> scale_factor = std::nullopt;
    std::optional<double> add_offset   = std::nullopt;
    void determine_scale_and_offset();

public:
    StorageDataType data_type = StorageDataType::Short;
    NetCDFFile()              = default;
    NetCDFFile( const Topography & topography, Output output ) : TopographyFile::TopographyFile( topography, output ){};
    // NetCDFFile( const std::filesystem::path & path, const std::optional<TopographyCrop> & crop = std::nullopt );
    void save( const std::filesystem::path & path ) override;
};

} // namespace Flowy
