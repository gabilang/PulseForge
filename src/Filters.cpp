#include "bsp/Filters.hpp"

#include <algorithm>
#include <cmath>

namespace bsp {

static constexpr double kPi = 3.14159265358979323846;

double Biquad::process(double x) {
    // Transposed direct form II.
    double y = b0 * x + z1;
    z1 = b1 * x - a1 * y + z2;
    z2 = b2 * x - a2 * y;
    return y;
}

// --- RBJ audio-EQ-cookbook biquads (normalised so a0 = 1) ---

static Biquad fromRaw(double b0, double b1, double b2, double a0, double a1, double a2) {
    Biquad q;
    q.b0 = b0 / a0; q.b1 = b1 / a0; q.b2 = b2 / a0;
    q.a1 = a1 / a0; q.a2 = a2 / a0;
    return q;
}

Biquad biquadLowpass(double fc, double fs, double Q) {
    double w0 = 2.0 * kPi * fc / fs;
    double c = std::cos(w0), s = std::sin(w0);
    double alpha = s / (2.0 * Q);
    double b1 = 1.0 - c;
    double b0 = b1 / 2.0, b2 = b1 / 2.0;
    return fromRaw(b0, b1, b2, 1.0 + alpha, -2.0 * c, 1.0 - alpha);
}

Biquad biquadHighpass(double fc, double fs, double Q) {
    double w0 = 2.0 * kPi * fc / fs;
    double c = std::cos(w0), s = std::sin(w0);
    double alpha = s / (2.0 * Q);
    double b0 = (1.0 + c) / 2.0, b1 = -(1.0 + c), b2 = (1.0 + c) / 2.0;
    return fromRaw(b0, b1, b2, 1.0 + alpha, -2.0 * c, 1.0 - alpha);
}

Biquad biquadBandpass(double fc, double fs, double Q) {
    double w0 = 2.0 * kPi * fc / fs;
    double c = std::cos(w0), s = std::sin(w0);
    double alpha = s / (2.0 * Q);
    double b0 = alpha, b1 = 0.0, b2 = -alpha;
    return fromRaw(b0, b1, b2, 1.0 + alpha, -2.0 * c, 1.0 - alpha);
}

Biquad biquadNotch(double fc, double fs, double Q) {
    double w0 = 2.0 * kPi * fc / fs;
    double c = std::cos(w0), s = std::sin(w0);
    double alpha = s / (2.0 * Q);
    double b0 = 1.0, b1 = -2.0 * c, b2 = 1.0;
    return fromRaw(b0, b1, b2, 1.0 + alpha, -2.0 * c, 1.0 - alpha);
}

std::vector<double> applyBiquad(Biquad bq, const std::vector<double>& x) {
    bq.reset();
    std::vector<double> y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y[i] = bq.process(x[i]);
    return y;
}

std::vector<double> filtfilt(const std::vector<double>& x, std::vector<Biquad> chain) {
    const std::size_t N = x.size();
    if (N == 0 || chain.empty()) return x;

    std::size_t pad = std::min<std::size_t>(N - 1, 3 * chain.size() * 10 + 12);
    std::vector<double> ext;
    ext.reserve(N + 2 * pad);
    for (std::size_t i = 0; i < pad; ++i) ext.push_back(2.0 * x[0] - x[pad - i]);     // odd reflection
    for (double v : x) ext.push_back(v);
    for (std::size_t i = 0; i < pad; ++i) ext.push_back(2.0 * x[N - 1] - x[N - 2 - i]);

    auto forward = [&](std::vector<double> in) {
        for (Biquad bq : chain) {        // copy: fresh state per pass
            bq.reset();
            for (double& s : in) s = bq.process(s);
        }
        return in;
    };

    std::vector<double> y = forward(ext);
    std::reverse(y.begin(), y.end());
    y = forward(y);
    std::reverse(y.begin(), y.end());

    return std::vector<double>(y.begin() + pad, y.begin() + pad + N);
}

// --- FIR ---

static double sinc(double x) {
    if (std::fabs(x) < 1e-12) return 1.0;
    double a = kPi * x;
    return std::sin(a) / a;
}

static int forceOdd(int n) { return (n % 2 == 0) ? n + 1 : n; }

std::vector<double> designLowpassFIR(int numTaps, double fc, double fs, WindowType wt) {
    numTaps = forceOdd(std::max(numTaps, 3));
    std::vector<double> w = makeWindow(wt, static_cast<std::size_t>(numTaps));
    std::vector<double> h(numTaps);
    double fcn = fc / fs;                         // normalised cutoff (cycles/sample)
    double mid = (numTaps - 1) / 2.0;
    double sum = 0.0;
    for (int n = 0; n < numTaps; ++n) {
        double m = n - mid;
        h[n] = 2.0 * fcn * sinc(2.0 * fcn * m) * w[n];
        sum += h[n];
    }
    for (double& v : h) v /= sum;                 // unity DC gain
    return h;
}

std::vector<double> designHighpassFIR(int numTaps, double fc, double fs, WindowType wt) {
    numTaps = forceOdd(std::max(numTaps, 3));
    std::vector<double> lp = designLowpassFIR(numTaps, fc, fs, wt);  // spectral inversion
    int mid = (numTaps - 1) / 2;
    for (int n = 0; n < numTaps; ++n) lp[n] = -lp[n];
    lp[mid] += 1.0;
    return lp;
}

std::vector<double> designBandpassFIR(int numTaps, double f1, double f2, double fs, WindowType wt) {
    numTaps = forceOdd(std::max(numTaps, 3));
    std::vector<double> low = designLowpassFIR(numTaps, f2, fs, wt);
    std::vector<double> highCut = designLowpassFIR(numTaps, f1, fs, wt);
    std::vector<double> h(numTaps);
    for (int n = 0; n < numTaps; ++n) h[n] = low[n] - highCut[n];   // LP(f2) - LP(f1)
    return h;
}

std::vector<double> applyFIR(const std::vector<double>& h, const std::vector<double>& x) {
    const std::size_t N = x.size(), M = h.size();
    if (N == 0 || M == 0) return x;
    const long delay = static_cast<long>(M - 1) / 2;   // symmetric FIR group delay
    std::vector<double> y(N, 0.0);
    for (std::size_t n = 0; n < N; ++n) {
        double acc = 0.0;
        for (std::size_t k = 0; k < M; ++k) {
            long idx = static_cast<long>(n) + delay - static_cast<long>(k);
            if (idx >= 0 && idx < static_cast<long>(N)) acc += h[k] * x[static_cast<std::size_t>(idx)];
        }
        y[n] = acc;
    }
    return y;
}

std::vector<double> movingAverage(const std::vector<double>& x, int window) {
    if (window < 1) window = 1;
    const std::size_t N = x.size();
    std::vector<double> y(N, 0.0);
    int half = window / 2;
    double sum = 0.0;
    // Naive centred average (clear over clever for this size).
    for (std::size_t n = 0; n < N; ++n) {
        long lo = static_cast<long>(n) - half;
        long hi = static_cast<long>(n) + half;
        if (lo < 0) lo = 0;
        if (hi >= static_cast<long>(N)) hi = static_cast<long>(N) - 1;
        sum = 0.0;
        for (long i = lo; i <= hi; ++i) sum += x[static_cast<std::size_t>(i)];
        y[n] = sum / static_cast<double>(hi - lo + 1);
    }
    return y;
}

std::vector<double> notchFilter(const std::vector<double>& x, double fs, double f0, double Q) {
    return filtfilt(x, {biquadNotch(f0, fs, Q)});
}

std::vector<double> bandpassZeroPhase(const std::vector<double>& x, double fs,
                                      double low, double high, double Q) {
    std::vector<Biquad> chain;
    if (low > 0.0)        chain.push_back(biquadHighpass(low, fs, Q));
    if (high < fs / 2.0)  chain.push_back(biquadLowpass(high, fs, Q));
    if (chain.empty()) return x;
    return filtfilt(x, chain);
}

} // namespace bsp
