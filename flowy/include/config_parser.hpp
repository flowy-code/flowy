#pragma once
// GPL v3 License
// Copyright 2023--present Flowy developers
#include "flowy/include/config.hpp"
#include <filesystem>

namespace Flowy::Config
{
InputParams parse_config( const std::filesystem::path & path );
void validate_settings( const InputParams & options );
} // namespace Flowy::Config
