#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/topography_file.hpp"
#include <optional>

namespace Flowy
{

class AscFile : public TopographyFile
{
private:
    std::string suffix = ".asc";

public:
    AscFile() = default;
    AscFile( const Topography & topography, OutputQuantitiy output )
            : TopographyFile::TopographyFile( topography, output ){};
    AscFile( const std::filesystem::path & path, const std::optional<TopographyCrop> & crop = std::nullopt );
    void save( const std::filesystem::path & path ) override;
};

} // namespace Flowy
