#include "bsp/viz/AsciiPlot.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <sstream>

namespace bsp {
namespace viz {

static constexpr int kLabelW = 10;  // left margin reserved for y-axis labels

std::string asciiPlot(const std::vector<double>& y, const AsciiOptions& opt) {
    std::ostringstream os;
    const int W = std::max(10, opt.width);
    const int H = std::max(3, opt.height);

    if (!opt.title.empty()) os << opt.title << "\n";
    if (y.empty()) { os << "(no data)\n"; return os.str(); }

    double dmin = *std::min_element(y.begin(), y.end());
    double dmax = *std::max_element(y.begin(), y.end());
    if (opt.baselineZero) { dmin = std::min(dmin, 0.0); dmax = std::max(dmax, 0.0); }
    if (dmax - dmin < 1e-12) dmax = dmin + 1.0;

    auto rowOf = [&](double v) {
        int r = static_cast<int>(std::lround((dmax - v) / (dmax - dmin) * (H - 1)));
        return std::min(std::max(r, 0), H - 1);
    };

    std::vector<std::string> grid(H, std::string(W, ' '));
    int zeroRow = (dmin <= 0.0 && dmax >= 0.0) ? rowOf(0.0) : -1;
    if (zeroRow >= 0) grid[zeroRow] = std::string(W, '-');

    const std::size_t N = y.size();
    for (int c = 0; c < W; ++c) {
        std::size_t a = static_cast<std::size_t>(static_cast<long long>(c) * N / W);
        std::size_t b = static_cast<std::size_t>(static_cast<long long>(c + 1) * N / W);
        if (b <= a) b = a + 1;
        if (b > N) b = N;
        double vmin = y[a], vmax = y[a];
        for (std::size_t i = a; i < b; ++i) { vmin = std::min(vmin, y[i]); vmax = std::max(vmax, y[i]); }

        int top, bot;
        if (opt.baselineZero) {
            top = std::min(rowOf(vmax), rowOf(0.0));
            bot = std::max(rowOf(vmax), rowOf(0.0));
        } else {
            top = rowOf(vmax);
            bot = rowOf(vmin);
        }
        for (int r = top; r <= bot; ++r) grid[r][c] = '*';
    }

    char buf[32];
    auto label = [&](double v) {
        std::snprintf(buf, sizeof buf, "%9.3g", v);
        return std::string(buf) + " ";
    };
    for (int r = 0; r < H; ++r) {
        std::string lab(kLabelW, ' ');
        if (r == 0) lab = label(dmax);
        else if (r == H - 1) lab = label(dmin);
        else if (r == zeroRow) lab = label(0.0);
        os << lab << "|" << grid[r] << "\n";
    }

    os << std::string(kLabelW, ' ') << "+" << std::string(W, '-') << "\n";

    if (opt.xMax > opt.xMin) {
        std::snprintf(buf, sizeof buf, "%.2f", opt.xMin);
        std::string left(buf);
        std::snprintf(buf, sizeof buf, "%.2f", opt.xMax);
        std::string right(buf);
        std::string line(kLabelW + 1 + W, ' ');
        for (std::size_t i = 0; i < left.size() && kLabelW + 1 + i < line.size(); ++i)
            line[kLabelW + 1 + i] = left[i];
        for (std::size_t i = 0; i < right.size(); ++i) {
            std::size_t pos = kLabelW + 1 + W - right.size() + i;
            if (pos < line.size()) line[pos] = right[i];
        }
        os << line << "\n";
        if (!opt.xlabel.empty()) {
            std::string l2(kLabelW + 1 + W, ' ');
            int start = kLabelW + 1 + (W - std::min<int>(opt.xlabel.size(), W)) / 2;
            for (std::size_t i = 0; i < opt.xlabel.size() && start + i < l2.size(); ++i)
                l2[start + i] = opt.xlabel[i];
            os << l2 << "\n";
        }
    }
    return os.str();
}

std::string asciiWaveform(const Signal& sig, int width, int height, const std::string& title) {
    AsciiOptions o;
    o.width = width;
    o.height = height;
    o.title = title.empty() ? (sig.label.empty() ? "waveform" : sig.label + " waveform") : title;
    o.xMin = 0.0;
    o.xMax = sig.duration();
    o.xlabel = "time (s)";
    return asciiPlot(sig.samples, o);
}

std::string asciiSpectrum(const PSD& psd, double fMax, int width, int height,
                          const std::string& title) {
    double maxF = fMax > 0.0 ? fMax : (psd.freqs.empty() ? 0.0 : psd.freqs.back());
    std::vector<double> y;
    for (std::size_t i = 0; i < psd.freqs.size(); ++i)
        if (psd.freqs[i] <= maxF) y.push_back(psd.power[i]);

    AsciiOptions o;
    o.width = width;
    o.height = height;
    o.baselineZero = true;
    o.title = title.empty() ? "power spectrum" : title;
    o.xMin = 0.0;
    o.xMax = maxF;
    o.xlabel = "frequency (Hz)";
    return asciiPlot(y, o);
}

} // namespace viz
} // namespace bsp
