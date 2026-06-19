#include "bsp/bio/EEG.hpp"

#include <algorithm>

#include "bsp/Statistics.hpp"

namespace bsp {
namespace bio {

EEGResult analyzeEEG(const Signal& sig, EEGBands bands) {
    EEGResult r;
    if (sig.empty()) return r;

    // Remove DC so very-low-frequency drift doesn't dominate.
    double m = stats::mean(sig.samples);
    std::vector<double> x(sig.samples.size());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = sig.samples[i] - m;

    // Welch PSD: ~2 s segments (power of two friendly) with 50% overlap.
    std::size_t seg = static_cast<std::size_t>(sig.fs * 2.0);
    if (seg < 64) seg = 64;
    r.psd = welch(x, sig.fs, seg, 0.5, WindowType::Hann);

    r.absolute.delta = bandPower(r.psd, bands.deltaLo, bands.deltaHi);
    r.absolute.theta = bandPower(r.psd, bands.thetaLo, bands.thetaHi);
    r.absolute.alpha = bandPower(r.psd, bands.alphaLo, bands.alphaHi);
    r.absolute.beta  = bandPower(r.psd, bands.betaLo,  bands.betaHi);
    r.absolute.gamma = bandPower(r.psd, bands.gammaLo, bands.gammaHi);
    r.absolute.total = r.absolute.delta + r.absolute.theta + r.absolute.alpha +
                       r.absolute.beta + r.absolute.gamma;

    double t = r.absolute.total > 0.0 ? r.absolute.total : 1.0;
    r.relative.delta = r.absolute.delta / t;
    r.relative.theta = r.absolute.theta / t;
    r.relative.alpha = r.absolute.alpha / t;
    r.relative.beta  = r.absolute.beta / t;
    r.relative.gamma = r.absolute.gamma / t;
    r.relative.total = 1.0;

    r.dominantFreq   = dominantFrequency(r.psd);
    r.spectralEdge95 = spectralEdgeFrequency(r.psd, 0.95);
    r.meanFreq       = meanFrequency(r.psd);
    return r;
}

} // namespace bio
} // namespace bsp
