#include "bsp/Statistics.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace bsp {
namespace stats {

double mean(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    double s = 0.0;
    for (double v : x) s += v;
    return s / static_cast<double>(x.size());
}

double variance(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    double m = mean(x), s = 0.0;
    for (double v : x) { double d = v - m; s += d * d; }
    return s / static_cast<double>(x.size() - 1);
}

double stddev(const std::vector<double>& x) { return std::sqrt(variance(x)); }

double rms(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    double s = 0.0;
    for (double v : x) s += v * v;
    return std::sqrt(s / static_cast<double>(x.size()));
}

double minValue(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    return *std::min_element(x.begin(), x.end());
}

double maxValue(const std::vector<double>& x) {
    if (x.empty()) return 0.0;
    return *std::max_element(x.begin(), x.end());
}

double median(std::vector<double> x) {
    if (x.empty()) return 0.0;
    std::sort(x.begin(), x.end());
    std::size_t n = x.size();
    return (n % 2) ? x[n / 2] : 0.5 * (x[n / 2 - 1] + x[n / 2]);
}

double skewness(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    double m = mean(x), sd = stddev(x);
    if (sd <= 0.0) return 0.0;
    double s = 0.0;
    for (double v : x) { double d = (v - m) / sd; s += d * d * d; }
    return s / static_cast<double>(x.size());
}

double kurtosis(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    double m = mean(x), sd = stddev(x);
    if (sd <= 0.0) return 0.0;
    double s = 0.0;
    for (double v : x) { double d = (v - m) / sd; s += d * d * d * d; }
    return s / static_cast<double>(x.size()) - 3.0;  // excess kurtosis
}

std::size_t zeroCrossings(const std::vector<double>& x) {
    std::size_t count = 0;
    for (std::size_t i = 1; i < x.size(); ++i) {
        if ((x[i - 1] < 0.0 && x[i] >= 0.0) || (x[i - 1] >= 0.0 && x[i] < 0.0)) ++count;
    }
    return count;
}

double zeroCrossingRate(const std::vector<double>& x) {
    if (x.size() < 2) return 0.0;
    return static_cast<double>(zeroCrossings(x)) / static_cast<double>(x.size() - 1);
}

} // namespace stats
} // namespace bsp
