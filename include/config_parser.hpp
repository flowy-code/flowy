#pragma once
#include "config.hpp"
#include <filesystem>

namespace Flowy::Config
{
InputParams parse_config( const std::filesystem::path & path );
void validate_settings( const InputParams & options );
} // namespace Flowy::Config