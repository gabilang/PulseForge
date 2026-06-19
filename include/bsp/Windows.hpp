#pragma once
#include <cstddef>
#include <string>
#include <vector>

namespace bsp {

enum class WindowType { Rectangular, Hann, Hamming, Blackman };

// Build a window of length N. N <= 1 returns a vector of ones.
std::vector<double> makeWindow(WindowType type, std::size_t N);

// Sum of squared window coefficients (used for PSD normalisation).
double windowPowerSum(const std::vector<double>& w);

WindowType windowFromString(const std::string& name);
const char* windowName(WindowType type);

} // namespace bsp
