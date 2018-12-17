#pragma once

#include "Config.hpp"

#include <optional>

namespace CommandLineParser
{

std::optional< Config > parse( int argc, char** argv );

} // namespace CommandLineParser
