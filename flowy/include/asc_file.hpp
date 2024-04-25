#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/topography_file.hpp"
#include <optional>
#include <string>

namespace Flowy
{

class AscFile : public TopographyFile
{
public:
    AscFile() = default;
    AscFile( const Topography & topography, OutputQuantitiy output )
            : TopographyFile::TopographyFile( topography, output )
    {
    }
    explicit AscFile( const std::filesystem::path & path, const std::optional<TopographyCrop> & crop = std::nullopt );

    std::string suffix() override
    {
        return ".asc";
    }

    void save( const std::filesystem::path & path ) override;
};

} // namespace Flowy
