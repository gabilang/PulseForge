#include "bsp/Generators.hpp"

#include <cmath>
#include <random>

namespace bsp {
namespace gen {

static constexpr double kPi = 3.14159265358979323846;

std::vector<double> sine(double freqHz, double amp, double fs, double durationSec, double phase) {
    std::size_t N = static_cast<std::size_t>(durationSec * fs);
    std::vector<double> y(N);
    for (std::size_t n = 0; n < N; ++n)
        y[n] = amp * std::sin(2.0 * kPi * freqHz * static_cast<double>(n) / fs + phase);
    return y;
}

std::vector<double> whiteNoise(double sigma, double fs, double durationSec, std::uint32_t seed) {
    std::size_t N = static_cast<std::size_t>(durationSec * fs);
    std::vector<double> y(N);
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, sigma);
    for (std::size_t n = 0; n < N; ++n) y[n] = nd(rng);
    return y;
}

static double gaussianBump(double t, double center, double amp, double width) {
    double d = (t - center) / width;
    return amp * std::exp(-0.5 * d * d);
}

Signal syntheticECG(double fs, double durationSec, double hr, double noise, std::uint32_t seed) {
    std::size_t N = static_cast<std::size_t>(durationSec * fs);
    std::vector<double> y(N, 0.0);
    double beatPeriod = 60.0 / hr;  // seconds between beats

    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, noise);

    for (std::size_t n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / fs;
        // Distance to the nearest beat centre.
        double phase = std::fmod(t, beatPeriod);
        double v = 0.0;
        // PQRST modelled as a sum of Gaussians relative to the R peak.
        v += gaussianBump(phase, beatPeriod * 0.0 - 0.0 + 0.20, 0.10, 0.025);  // P (wraps near prev beat)
        v += gaussianBump(phase, 0.36, -0.12, 0.012);                          // Q
        v += gaussianBump(phase, 0.40,  1.00, 0.011);                          // R
        v += gaussianBump(phase, 0.44, -0.22, 0.012);                          // S
        v += gaussianBump(phase, 0.58,  0.30, 0.040);                          // T
        y[n] = v + nd(rng);
    }
    return Signal(std::move(y), fs, "ECG");
}

Signal syntheticEEG(double fs, double durationSec, double noise, std::uint32_t seed) {
    std::size_t N = static_cast<std::size_t>(durationSec * fs);
    std::vector<double> y(N, 0.0);
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, noise);

    struct Comp { double f, a; };
    // Alpha (10 Hz) dominant, as in relaxed eyes-closed EEG.
    const Comp comps[] = {{2.0, 8.0}, {6.0, 6.0}, {10.0, 20.0}, {20.0, 5.0}, {40.0, 2.0}};
    for (std::size_t n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / fs;
        double v = 0.0;
        for (const auto& c : comps) v += c.a * std::sin(2.0 * kPi * c.f * t);
        y[n] = v + nd(rng);
    }
    return Signal(std::move(y), fs, "EEG");
}

Signal syntheticEMG(double fs, double durationSec, std::uint32_t seed) {
    std::size_t N = static_cast<std::size_t>(durationSec * fs);
    std::vector<double> y(N, 0.0);
    std::mt19937 rng(seed);
    std::normal_distribution<double> baseline(0.0, 0.02);
    std::normal_distribution<double> burst(0.0, 1.0);

    // Activation windows (start, end) in seconds.
    const double bursts[][2] = {{1.0, 2.0}, {4.0, 5.0}, {7.0, 8.5}};
    for (std::size_t n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / fs;
        double v = baseline(rng);
        for (const auto& b : bursts) {
            if (t >= b[0] && t <= b[1]) {
                // Smooth ramp envelope inside the burst.
                double mid = 0.5 * (b[0] + b[1]);
                double halfw = 0.5 * (b[1] - b[0]);
                double env = std::max(0.0, 1.0 - std::fabs(t - mid) / halfw);
                v += env * burst(rng);
            }
        }
        y[n] = v;
    }
    return Signal(std::move(y), fs, "EMG");
}

Signal syntheticEOG(double fs, double durationSec, double blinkRate, std::uint32_t seed) {
    std::size_t N = static_cast<std::size_t>(durationSec * fs);
    std::vector<double> y(N, 0.0);
    std::mt19937 rng(seed);
    std::normal_distribution<double> nd(0.0, 0.03);

    double blinkPeriod = 60.0 / blinkRate;  // seconds between blinks
    for (std::size_t n = 0; n < N; ++n) {
        double t = static_cast<double>(n) / fs;
        // Slow baseline drift.
        double v = 0.15 * std::sin(2.0 * kPi * 0.05 * t);
        // Periodic blink bumps.
        double phase = std::fmod(t, blinkPeriod);
        v += gaussianBump(phase, blinkPeriod * 0.5, 1.0, 0.06);
        y[n] = v + nd(rng);
    }
    return Signal(std::move(y), fs, "EOG");
}

} // namespace gen
} // namespace bsp
