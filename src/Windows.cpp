#include "bsp/Windows.hpp"

#include <algorithm>
#include <cmath>

namespace bsp {

static constexpr double kPi = 3.14159265358979323846;

std::vector<double> makeWindow(WindowType type, std::size_t N) {
    std::vector<double> w(N, 1.0);
    if (N <= 1) return w;
    const double M = static_cast<double>(N - 1);
    for (std::size_t n = 0; n < N; ++n) {
        const double r = static_cast<double>(n) / M;  // 0..1
        switch (type) {
            case WindowType::Rectangular:
                w[n] = 1.0;
                break;
            case WindowType::Hann:
                w[n] = 0.5 - 0.5 * std::cos(2.0 * kPi * r);
                break;
            case WindowType::Hamming:
                w[n] = 0.54 - 0.46 * std::cos(2.0 * kPi * r);
                break;
            case WindowType::Blackman:
                w[n] = 0.42 - 0.5 * std::cos(2.0 * kPi * r) + 0.08 * std::cos(4.0 * kPi * r);
                break;
        }
    }
    return w;
}

double windowPowerSum(const std::vector<double>& w) {
    double s = 0.0;
    for (double v : w) s += v * v;
    return s;
}

WindowType windowFromString(const std::string& name) {
    std::string n;
    n.reserve(name.size());
    for (char c : name) n.push_back(static_cast<char>(std::tolower(c)));
    if (n == "rect" || n == "rectangular" || n == "none") return WindowType::Rectangular;
    if (n == "hamming") return WindowType::Hamming;
    if (n == "blackman") return WindowType::Blackman;
    return WindowType::Hann;
}

const char* windowName(WindowType type) {
    switch (type) {
        case WindowType::Rectangular: return "rectangular";
        case WindowType::Hann:        return "hann";
        case WindowType::Hamming:     return "hamming";
        case WindowType::Blackman:    return "blackman";
    }
    return "hann";
}

} // namespace bsp
