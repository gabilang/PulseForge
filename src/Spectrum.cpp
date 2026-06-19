#include "bsp/Spectrum.hpp"

#include <algorithm>
#include <cmath>

#include "bsp/FFT.hpp"

namespace bsp {

// Modified periodogram of one windowed segment, one-sided, scaled to a
// proper density (units^2/Hz).
static void accumulateSegment(const double* seg, std::size_t segLen, double fs,
                              const std::vector<double>& win, double winPow,
                              std::vector<double>& acc) {
    const std::size_t nfft = nextPowerOfTwo(segLen);
    CArray a(nfft, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < segLen; ++i) a[i] = Complex(seg[i] * win[i], 0.0);
    fft(a);

    const std::size_t half = nfft / 2 + 1;
    const double scale = 1.0 / (fs * winPow);
    for (std::size_t k = 0; k < half; ++k) {
        double p = std::norm(a[k]) * scale;          // |X|^2 / (fs * sum(w^2))
        if (k != 0 && k != nfft / 2) p *= 2.0;        // one-sided (fold negatives)
        acc[k] += p;
    }
}

PSD welch(const std::vector<double>& x, double fs, std::size_t segLen,
          double overlapFrac, WindowType winType) {
    PSD out;
    if (x.empty() || fs <= 0.0) return out;
    if (segLen == 0 || segLen > x.size()) segLen = x.size();
    if (overlapFrac < 0.0) overlapFrac = 0.0;
    if (overlapFrac >= 1.0) overlapFrac = 0.9;

    const std::size_t nfft = nextPowerOfTwo(segLen);
    const std::size_t half = nfft / 2 + 1;
    const std::vector<double> win = makeWindow(winType, segLen);
    const double winPow = windowPowerSum(win);

    std::size_t step = static_cast<std::size_t>(segLen * (1.0 - overlapFrac));
    if (step == 0) step = 1;

    std::vector<double> acc(half, 0.0);
    std::size_t nSeg = 0;
    for (std::size_t start = 0; start + segLen <= x.size(); start += step) {
        accumulateSegment(&x[start], segLen, fs, win, winPow, acc);
        ++nSeg;
    }
    if (nSeg == 0) {  // signal shorter than one full segment
        std::vector<double> seg = x;
        seg.resize(segLen, 0.0);
        accumulateSegment(seg.data(), segLen, fs, win, winPow, acc);
        nSeg = 1;
    }

    out.resolution = fs / static_cast<double>(nfft);
    out.freqs.resize(half);
    out.power.resize(half);
    for (std::size_t k = 0; k < half; ++k) {
        out.freqs[k] = static_cast<double>(k) * out.resolution;
        out.power[k] = acc[k] / static_cast<double>(nSeg);
    }
    return out;
}

PSD periodogram(const std::vector<double>& x, double fs, WindowType win) {
    return welch(x, fs, x.size(), 0.0, win);
}

double bandPower(const PSD& psd, double f1, double f2) {
    if (psd.freqs.size() < 2) return 0.0;
    if (f1 > f2) std::swap(f1, f2);
    double total = 0.0;
    for (std::size_t k = 1; k < psd.freqs.size(); ++k) {
        double fa = psd.freqs[k - 1], fb = psd.freqs[k];
        if (fb <= f1 || fa >= f2) continue;          // bin fully outside band
        double lo = std::max(fa, f1), hi = std::min(fb, f2);
        // Linear interpolation of power density at the clipped edges.
        double pa = psd.power[k - 1], pb = psd.power[k];
        double slope = (pb - pa) / (fb - fa);
        double plo = pa + slope * (lo - fa);
        double phi = pa + slope * (hi - fa);
        total += 0.5 * (plo + phi) * (hi - lo);       // trapezoid over [lo,hi]
    }
    return total;
}

double totalPower(const PSD& psd) {
    if (psd.freqs.empty()) return 0.0;
    return bandPower(psd, psd.freqs.front(), psd.freqs.back());
}

double dominantFrequency(const PSD& psd) {
    if (psd.power.empty()) return 0.0;
    std::size_t best = 0;
    double bestP = -1.0;
    for (std::size_t k = 0; k < psd.power.size(); ++k) {
        if (psd.power[k] > bestP) { bestP = psd.power[k]; best = k; }
    }
    return psd.freqs[best];
}

double spectralEdgeFrequency(const PSD& psd, double edgeFraction) {
    double total = totalPower(psd);
    if (total <= 0.0) return 0.0;
    double target = total * edgeFraction;
    double cum = 0.0;
    for (std::size_t k = 1; k < psd.freqs.size(); ++k) {
        double seg = 0.5 * (psd.power[k - 1] + psd.power[k]) * (psd.freqs[k] - psd.freqs[k - 1]);
        if (cum + seg >= target) {
            double need = target - cum;
            double frac = seg > 0.0 ? need / seg : 0.0;
            return psd.freqs[k - 1] + frac * (psd.freqs[k] - psd.freqs[k - 1]);
        }
        cum += seg;
    }
    return psd.freqs.back();
}

double meanFrequency(const PSD& psd) {
    double num = 0.0, den = 0.0;
    for (std::size_t k = 0; k < psd.freqs.size(); ++k) {
        num += psd.freqs[k] * psd.power[k];
        den += psd.power[k];
    }
    return den > 0.0 ? num / den : 0.0;
}

double medianFrequency(const PSD& psd) {
    return spectralEdgeFrequency(psd, 0.5);
}

} // namespace bsp
