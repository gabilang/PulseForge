#include "bsp/bio/ECG.hpp"

#include <algorithm>
#include <cmath>

#include "bsp/Filters.hpp"
#include "bsp/Statistics.hpp"

namespace bsp {
namespace bio {

ECGResult analyzeECG(const Signal& sig) {
    ECGResult r;
    if (sig.size() < static_cast<std::size_t>(sig.fs)) return r;  // need ~1 s
    const double fs = sig.fs;

    // 1. Band-pass 5-15 Hz to isolate QRS energy (zero-phase).
    std::vector<double> bp = bandpassZeroPhase(sig.samples, fs, 5.0, 15.0);

    // 2. Five-point derivative emphasises the steep QRS slope.
    std::vector<double> deriv(bp.size(), 0.0);
    for (std::size_t n = 4; n < bp.size(); ++n)
        deriv[n] = (2.0 * bp[n] + bp[n - 1] - bp[n - 3] - 2.0 * bp[n - 4]) / 8.0;

    // 3. Square.
    std::vector<double> sq(deriv.size());
    for (std::size_t n = 0; n < deriv.size(); ++n) sq[n] = deriv[n] * deriv[n];

    // 4. Moving-window integration (~150 ms).
    int win = std::max(1, static_cast<int>(0.150 * fs));
    std::vector<double> mwi = movingAverage(sq, win);

    // 5. Threshold + refractory-period peak picking on the integrated signal.
    double peak = stats::maxValue(mwi);
    if (peak <= 0.0) return r;
    double thr = 0.35 * peak;
    std::size_t refractory = static_cast<std::size_t>(0.25 * fs);
    std::size_t searchBack = static_cast<std::size_t>(0.20 * fs);

    std::vector<std::size_t> peaks;
    std::size_t i = 0;
    while (i < mwi.size()) {
        if (mwi[i] > thr) {
            // Find the local maximum within the supra-threshold region.
            std::size_t regionMax = i;
            while (i < mwi.size() && mwi[i] > thr) {
                if (mwi[i] > mwi[regionMax]) regionMax = i;
                ++i;
            }
            // Refine to the true R peak on the band-passed trace just before
            // the integrator peak (the integrator introduces a lag).
            std::size_t lo = regionMax > searchBack ? regionMax - searchBack : 0;
            std::size_t hi = std::min(regionMax + 1, bp.size());
            std::size_t rpos = lo;
            for (std::size_t k = lo; k < hi; ++k)
                if (std::fabs(bp[k]) > std::fabs(bp[rpos])) rpos = k;

            if (peaks.empty() || rpos - peaks.back() >= refractory)
                peaks.push_back(rpos);
            else if (std::fabs(bp[rpos]) > std::fabs(bp[peaks.back()]))
                peaks.back() = rpos;  // keep the stronger of two too-close peaks
        } else {
            ++i;
        }
    }

    r.rPeaks = peaks;
    if (peaks.size() < 2) return r;

    // R-R intervals and heart-rate / HRV metrics.
    std::vector<double> rrSec, hrInst;
    for (std::size_t k = 1; k < peaks.size(); ++k) {
        double rr = static_cast<double>(peaks[k] - peaks[k - 1]) / fs;
        rrSec.push_back(rr);
        hrInst.push_back(60.0 / rr);
    }
    r.rrIntervals = rrSec;
    r.meanHR = stats::mean(hrInst);
    r.minHR = stats::minValue(hrInst);
    r.maxHR = stats::maxValue(hrInst);

    // SDNN (ms).
    std::vector<double> rrMs;
    for (double s : rrSec) rrMs.push_back(s * 1000.0);
    r.sdnn = stats::stddev(rrMs);

    // RMSSD and pNN50.
    if (rrMs.size() >= 2) {
        double sumSq = 0.0;
        int nn50 = 0;
        for (std::size_t k = 1; k < rrMs.size(); ++k) {
            double d = rrMs[k] - rrMs[k - 1];
            sumSq += d * d;
            if (std::fabs(d) > 50.0) ++nn50;
        }
        r.rmssd = std::sqrt(sumSq / static_cast<double>(rrMs.size() - 1));
        r.pnn50 = 100.0 * static_cast<double>(nn50) / static_cast<double>(rrMs.size() - 1);
    }
    return r;
}

} // namespace bio
} // namespace bsp
