#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/definitions.hpp"
#include "flowy/include/topography_file.hpp"
#include <optional>

namespace Flowy
{

class AscFile : public TopographyFile
{
public:
    AscFile() = default;
    AscFile( const std::filesystem::path & path, std::optional<AscCrop> crop = std::nullopt );
    void save( const std::filesystem::path & path ) override;
};

} // namespace Flowy
