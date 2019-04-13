#pragma once

#include "Config.hpp"

#include <optional>

namespace GanonClassify::CommandLineParser
{

std::optional< Config > parse( int argc, char** argv );

} // namespace GanonClassify::CommandLineParser
