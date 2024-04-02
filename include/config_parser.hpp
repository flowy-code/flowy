#pragma once
#include "config.hpp"
#include <filesystem>

namespace Flowtastic::Config
{
InputParams parse_config( const std::filesystem::path & path );
}