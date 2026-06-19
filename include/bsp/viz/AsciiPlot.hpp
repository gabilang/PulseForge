#pragma once
#include <string>
#include <vector>

#include "bsp/Signal.hpp"
#include "bsp/Spectrum.hpp"

namespace bsp {
namespace viz {

// Render plots as monospace text suitable for printing to a terminal.
struct AsciiOptions {
    int width = 78;          // plot area width in characters
    int height = 16;         // plot area height in rows
    std::string title;
    std::string xlabel;
    double xMin = 0.0;
    double xMax = 0.0;       // if xMax > xMin, draw an x-axis range
    bool baselineZero = false; // draw bars from a zero baseline (for spectra)
};

std::string asciiPlot(const std::vector<double>& y, const AsciiOptions& opt = {});

std::string asciiWaveform(const Signal& sig, int width = 78, int height = 16,
                          const std::string& title = "");

std::string asciiSpectrum(const PSD& psd, double fMax, int width = 78, int height = 16,
                          const std::string& title = "");

} // namespace viz
} // namespace bsp
