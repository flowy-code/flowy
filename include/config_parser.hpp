#pragma once
#include "config.hpp"
#include <filesystem>

namespace Flowy::Config
{
InputParams parse_config( const std::filesystem::path & path );
}