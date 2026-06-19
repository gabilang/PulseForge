#include "bsp/bio/EOG.hpp"

#include <algorithm>
#include <cmath>

#include "bsp/Filters.hpp"
#include "bsp/Statistics.hpp"

namespace bsp {
namespace bio {

EOGResult analyzeEOG(const Signal& sig) {
    EOGResult r;
    if (sig.empty()) return r;
    const double fs = sig.fs;

    // Low-pass at 10 Hz: blinks and saccades are slow, large deflections.
    double cutoff = std::min(10.0, fs * 0.4);
    std::vector<double> lp = bandpassZeroPhase(sig.samples, fs, 0.1, cutoff);

    // Threshold relative to signal spread.
    double m = stats::mean(lp);
    double sd = stats::stddev(lp);
    double thr = m + 2.5 * sd;

    std::size_t refractory = static_cast<std::size_t>(0.20 * fs);  // blinks ≥200 ms apart
    std::vector<std::size_t> blinks;
    double ampSum = 0.0;

    std::size_t i = 0;
    while (i < lp.size()) {
        if (lp[i] > thr) {
            std::size_t regionMax = i;
            while (i < lp.size() && lp[i] > thr) {
                if (lp[i] > lp[regionMax]) regionMax = i;
                ++i;
            }
            if (blinks.empty() || regionMax - blinks.back() >= refractory) {
                blinks.push_back(regionMax);
                ampSum += lp[regionMax];
            }
        } else {
            ++i;
        }
    }

    r.blinks = blinks;
    r.blinkCount = static_cast<int>(blinks.size());
    double durMin = sig.duration() / 60.0;
    r.blinkRatePerMin = durMin > 0.0 ? r.blinkCount / durMin : 0.0;
    r.meanBlinkAmplitude = blinks.empty() ? 0.0 : ampSum / static_cast<double>(blinks.size());
    return r;
}

} // namespace bio
} // namespace bsp
