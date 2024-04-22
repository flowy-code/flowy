#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers

#include "flowy/include/topography_file.hpp"

namespace Flowy
{

class NetCDFFile : public TopographyFile
{
public:
    NetCDFFile() = default;
    NetCDFFile( const Topography & topography, Output output ) : TopographyFile::TopographyFile( topography, output ){};
    // NetCDFFile( const std::filesystem::path & path, const std::optional<TopographyCrop> & crop = std::nullopt );
    void save( const std::filesystem::path & path ) override;
};

} // namespace Flowy
