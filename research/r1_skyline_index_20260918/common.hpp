#pragma once
// Stage-2 shared compilation marker.  The canonical-tree implementation is
// intentionally retained verbatim in count.cpp and is included by index.cpp
// (with its entry point renamed), so both executables compile the identical
// stage-1 construction.
#include "../r1_terminal_20260918/terminal.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
