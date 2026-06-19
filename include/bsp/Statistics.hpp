#pragma once
#include <cstddef>
#include <vector>

namespace bsp {
namespace stats {

double mean(const std::vector<double>& x);
double variance(const std::vector<double>& x);   // sample variance (N-1)
double stddev(const std::vector<double>& x);
double rms(const std::vector<double>& x);
double minValue(const std::vector<double>& x);
double maxValue(const std::vector<double>& x);
double median(std::vector<double> x);             // by value (sorts a copy)
double skewness(const std::vector<double>& x);
double kurtosis(const std::vector<double>& x);    // excess kurtosis
std::size_t zeroCrossings(const std::vector<double>& x);
double zeroCrossingRate(const std::vector<double>& x);

} // namespace stats
} // namespace bsp
