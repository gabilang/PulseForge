#pragma once
#include <cstddef>
#include <vector>

#include "bsp/Windows.hpp"

namespace bsp {

// One-sided power spectral density estimate.
struct PSD {
    std::vector<double> freqs;  // frequency bins (Hz)
    std::vector<double> power;  // power density at each bin (units^2/Hz)
    double resolution = 0.0;    // frequency resolution (Hz/bin)
};

// Welch's method: average modified periodograms over overlapping segments.
// segLen is the segment length in samples; overlapFrac in [0,1).
PSD welch(const std::vector<double>& x, double fs, std::size_t segLen,
          double overlapFrac = 0.5, WindowType win = WindowType::Hann);

// Single-segment periodogram over the whole signal.
PSD periodogram(const std::vector<double>& x, double fs, WindowType win = WindowType::Hann);

// Integrate the PSD over [f1, f2] using the trapezoidal rule.
double bandPower(const PSD& psd, double f1, double f2);

// Total power across the whole spectrum.
double totalPower(const PSD& psd);

// Frequency of the largest spectral peak.
double dominantFrequency(const PSD& psd);

// Spectral edge frequency: the frequency below which `edgeFraction`
// (e.g. 0.95) of the total power lies.
double spectralEdgeFrequency(const PSD& psd, double edgeFraction);

// Power-weighted mean frequency.
double meanFrequency(const PSD& psd);

// Median (equal-power) frequency.
double medianFrequency(const PSD& psd);

} // namespace bsp
