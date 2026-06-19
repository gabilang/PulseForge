#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include "bsp/Signal.hpp"
#include "bsp/Spectrum.hpp"

namespace bsp {
namespace viz {

// Render plots as standalone SVG documents (open in any browser/viewer).
struct SvgOptions {
    int width = 900;
    int height = 300;
    std::string title;
    std::string xlabel = "time (s)";
    std::string lineColor = "#2563eb";
    std::string markerColor = "#dc2626";
    bool baselineZero = false;  // anchor the y-range at zero
    bool fillArea = false;      // shade under the curve
};

// Generic series plot: y sampled every dt x-units. `markers` are sample indices
// drawn as vertical reference lines + dots (e.g. ECG R-peaks, EOG blinks).
std::string svgSeries(const std::vector<double>& y, double dt, const SvgOptions& opt,
                      const std::vector<std::size_t>& markers = {});

std::string svgWaveform(const Signal& sig, const SvgOptions& opt = {},
                        const std::vector<std::size_t>& markers = {});

std::string svgSpectrum(const PSD& psd, double fMax, const SvgOptions& opt = {});

void writeSvg(const std::string& path, const std::string& content);

} // namespace viz
} // namespace bsp
