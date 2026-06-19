#pragma once
#include <string>
#include <vector>

#include "bsp/Signal.hpp"

namespace bsp {
namespace io {

// Load a single numeric column (0-based) from a CSV/whitespace file.
// Lines that fail to parse the requested column are skipped; a leading
// non-numeric header row is tolerated.
std::vector<double> loadColumn(const std::string& path, std::size_t column = 0);

// Load a column as a Signal with the given sample rate.
Signal loadSignal(const std::string& path, double fs, std::size_t column = 0,
                  const std::string& label = "");

// Write a signal as two columns: time,value.
void saveSignal(const std::string& path, const Signal& sig);

// Write arbitrary named columns (all assumed equal length of the first).
void saveColumns(const std::string& path, const std::vector<std::string>& headers,
                 const std::vector<std::vector<double>>& columns);

} // namespace io
} // namespace bsp
