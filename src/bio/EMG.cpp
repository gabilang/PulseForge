#include "bsp/bio/EMG.hpp"

#include <algorithm>
#include <cmath>

#include "bsp/Filters.hpp"
#include "bsp/Spectrum.hpp"
#include "bsp/Statistics.hpp"

namespace bsp {
namespace bio {

EMGResult analyzeEMG(const Signal& sig) {
    EMGResult r;
    if (sig.empty()) return r;
    const double fs = sig.fs;

    // Band-pass to the EMG band (20 Hz .. min(450, Nyquist*0.9)).
    double high = std::min(450.0, fs * 0.45);
    std::vector<double> x = bandpassZeroPhase(sig.samples, fs, 20.0, high);

    // RMS envelope: sliding root-mean-square over ~100 ms.
    int win = std::max(1, static_cast<int>(0.100 * fs));
    std::vector<double> sq(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) sq[i] = x[i] * x[i];
    std::vector<double> ms = movingAverage(sq, win);
    r.envelope.resize(ms.size());
    for (std::size_t i = 0; i < ms.size(); ++i) r.envelope[i] = std::sqrt(ms[i]);

    r.rms = stats::rms(x);
    r.maxEnvelope = stats::maxValue(r.envelope);

    // Onset detection: threshold at baseline + k*sigma, estimated from the
    // quietest 20% of the envelope.
    std::vector<double> sorted = r.envelope;
    std::sort(sorted.begin(), sorted.end());
    std::size_t q = std::max<std::size_t>(1, sorted.size() / 5);
    std::vector<double> quiet(sorted.begin(), sorted.begin() + q);
    double base = stats::mean(quiet);
    double baseSd = stats::stddev(quiet);
    double thr = base + 3.0 * baseSd + 0.10 * (r.maxEnvelope - base);

    std::size_t minLen = static_cast<std::size_t>(0.05 * fs);  // ignore <50 ms blips
    bool active = false;
    std::size_t startIdx = 0;
    for (std::size_t i = 0; i < r.envelope.size(); ++i) {
        if (!active && r.envelope[i] > thr) { active = true; startIdx = i; }
        else if (active && r.envelope[i] <= thr) {
            if (i - startIdx >= minLen)
                r.activations.push_back({startIdx / fs, i / fs});
            active = false;
        }
    }
    if (active && r.envelope.size() - startIdx >= minLen)
        r.activations.push_back({startIdx / fs, (r.envelope.size() - 1) / fs});

    // Spectral fatigue indices.
    PSD psd = welch(x, fs, std::min<std::size_t>(x.size(), static_cast<std::size_t>(fs)), 0.5,
                    WindowType::Hann);
    r.meanFrequency = meanFrequency(psd);
    r.medianFrequency = medianFrequency(psd);
    return r;
}

} // namespace bio
} // namespace bsp
