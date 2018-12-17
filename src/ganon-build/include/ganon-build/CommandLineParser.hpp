#pragma once

#include "Config.hpp"

#include <optional>

namespace GanonBuild::CommandLineParser
{

std::optional< Config > parse( int argc, char** argv );

} // namespace GanonBuild::CommandLineParser
