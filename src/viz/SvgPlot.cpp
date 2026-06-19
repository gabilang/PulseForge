#include "bsp/viz/SvgPlot.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace bsp {
namespace viz {

static std::string num(double v) {
    std::ostringstream o;
    o << v;
    return o.str();
}

std::string svgSeries(const std::vector<double>& y, double dt, const SvgOptions& opt,
                      const std::vector<std::size_t>& markers) {
    const int W = opt.width, H = opt.height;
    const int mL = 60, mR = 20, mB = 40;
    const int mT = opt.title.empty() ? 20 : 40;
    const double pw = W - mL - mR;
    const double ph = H - mT - mB;

    std::ostringstream s;
    s << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H
      << "\" viewBox=\"0 0 " << W << " " << H << "\" font-family=\"sans-serif\" font-size=\"12\">\n";
    s << "<rect x=\"0\" y=\"0\" width=\"" << W << "\" height=\"" << H << "\" fill=\"#ffffff\"/>\n";
    if (!opt.title.empty())
        s << "<text x=\"" << W / 2 << "\" y=\"22\" text-anchor=\"middle\" font-size=\"15\" "
             "font-weight=\"bold\">" << opt.title << "</text>\n";

    if (y.empty()) {
        s << "<text x=\"" << W / 2 << "\" y=\"" << H / 2 << "\" text-anchor=\"middle\">no data</text>\n</svg>\n";
        return s.str();
    }

    double dmin = *std::min_element(y.begin(), y.end());
    double dmax = *std::max_element(y.begin(), y.end());
    if (opt.baselineZero) { dmin = std::min(dmin, 0.0); dmax = std::max(dmax, 0.0); }
    if (dmax - dmin < 1e-12) dmax = dmin + 1.0;
    double xMax = dt * static_cast<double>(y.size() - 1);
    if (xMax <= 0.0) xMax = 1.0;

    auto X = [&](double t) { return mL + pw * (t / xMax); };
    auto Y = [&](double v) { return mT + ph * (dmax - v) / (dmax - dmin); };

    s << "<rect x=\"" << mL << "\" y=\"" << mT << "\" width=\"" << pw << "\" height=\"" << ph
      << "\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>\n";

    auto hline = [&](double v, const char* col, const char* dash) {
        double yy = Y(v);
        s << "<line x1=\"" << mL << "\" y1=\"" << yy << "\" x2=\"" << (mL + pw) << "\" y2=\"" << yy
          << "\" stroke=\"" << col << "\"";
        if (dash) s << " stroke-dasharray=\"" << dash << "\"";
        s << "/>\n";
    };
    hline(dmax, "#e2e8f0", nullptr);
    hline(dmin, "#e2e8f0", nullptr);
    if (dmin < 0.0 && dmax > 0.0) hline(0.0, "#94a3b8", "4 3");

    s << "<text x=\"" << (mL - 6) << "\" y=\"" << (Y(dmax) + 4) << "\" text-anchor=\"end\">"
      << num(dmax) << "</text>\n";
    s << "<text x=\"" << (mL - 6) << "\" y=\"" << (Y(dmin) + 4) << "\" text-anchor=\"end\">"
      << num(dmin) << "</text>\n";

    // Build the path, using min/max bucket decimation for long signals so peaks survive.
    std::ostringstream path;
    const std::size_t N = y.size();
    const std::size_t maxPts = 1500;
    if (N <= maxPts) {
        for (std::size_t i = 0; i < N; ++i)
            path << (i == 0 ? "M" : "L") << X(dt * i) << " " << Y(y[i]) << " ";
    } else {
        std::size_t buckets = maxPts / 2;
        bool first = true;
        for (std::size_t bk = 0; bk < buckets; ++bk) {
            std::size_t a = bk * N / buckets, e = (bk + 1) * N / buckets;
            if (e <= a) e = a + 1;
            if (e > N) e = N;
            std::size_t imin = a, imax = a;
            for (std::size_t i = a; i < e; ++i) {
                if (y[i] < y[imin]) imin = i;
                if (y[i] > y[imax]) imax = i;
            }
            std::size_t f1 = std::min(imin, imax), f2 = std::max(imin, imax);
            path << (first ? "M" : "L") << X(dt * f1) << " " << Y(y[f1]) << " "
                 << "L" << X(dt * f2) << " " << Y(y[f2]) << " ";
            first = false;
        }
    }

    if (opt.fillArea) {
        double baseY = Y(std::max(dmin, 0.0));
        s << "<path d=\"" << path.str() << "L" << X(xMax) << " " << baseY << " L" << X(0.0) << " "
          << baseY << " Z\" fill=\"" << opt.lineColor << "\" fill-opacity=\"0.15\" stroke=\"none\"/>\n";
    }
    s << "<path d=\"" << path.str() << "\" fill=\"none\" stroke=\"" << opt.lineColor
      << "\" stroke-width=\"1.2\"/>\n";

    for (std::size_t m : markers) {
        if (m >= N) continue;
        double xx = X(dt * m);
        s << "<line x1=\"" << xx << "\" y1=\"" << mT << "\" x2=\"" << xx << "\" y2=\"" << (mT + ph)
          << "\" stroke=\"" << opt.markerColor << "\" stroke-opacity=\"0.55\"/>\n";
        s << "<circle cx=\"" << xx << "\" cy=\"" << Y(y[m]) << "\" r=\"3\" fill=\"" << opt.markerColor
          << "\"/>\n";
    }

    s << "<text x=\"" << mL << "\" y=\"" << (mT + ph + 18) << "\" text-anchor=\"start\">0</text>\n";
    s << "<text x=\"" << (mL + pw) << "\" y=\"" << (mT + ph + 18) << "\" text-anchor=\"end\">"
      << num(xMax) << "</text>\n";
    if (!opt.xlabel.empty())
        s << "<text x=\"" << (mL + pw / 2) << "\" y=\"" << (H - 6) << "\" text-anchor=\"middle\">"
          << opt.xlabel << "</text>\n";

    // PulseForge brand mark: a small pulse glyph + wordmark, bottom-right.
    {
        double bx = mL + pw - 86, by = H - 8;
        s << "<g opacity=\"0.6\">\n";
        s << "<path d=\"M" << bx << " " << by << " l4 0 l2 -7 l3 11 l2 -4 l4 0\" fill=\"none\" stroke=\""
          << opt.markerColor << "\" stroke-width=\"1.4\"/>\n";
        s << "<text x=\"" << (bx + 22) << "\" y=\"" << (by + 2)
          << "\" font-size=\"11\" font-weight=\"bold\" fill=\"#475569\">PulseForge</text>\n";
        s << "</g>\n";
    }

    s << "</svg>\n";
    return s.str();
}

std::string svgWaveform(const Signal& sig, const SvgOptions& opt,
                        const std::vector<std::size_t>& markers) {
    SvgOptions o = opt;
    if (o.title.empty()) o.title = (sig.label.empty() ? "signal" : sig.label) + " waveform";
    if (o.xlabel.empty()) o.xlabel = "time (s)";
    double dt = sig.fs > 0.0 ? 1.0 / sig.fs : 1.0;
    return svgSeries(sig.samples, dt, o, markers);
}

std::string svgSpectrum(const PSD& psd, double fMax, const SvgOptions& opt) {
    SvgOptions o = opt;
    o.baselineZero = true;
    o.fillArea = true;
    if (o.title.empty()) o.title = "power spectrum";
    o.xlabel = "frequency (Hz)";

    double maxF = fMax > 0.0 ? fMax : (psd.freqs.empty() ? 1.0 : psd.freqs.back());
    std::vector<double> y;
    for (std::size_t i = 0; i < psd.freqs.size(); ++i)
        if (psd.freqs[i] <= maxF) y.push_back(psd.power[i]);
    double df = psd.resolution > 0.0 ? psd.resolution : 1.0;
    return svgSeries(y, df, o, {});
}

void writeSvg(const std::string& path, const std::string& content) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("Cannot write file: " + path);
    out << content;
}

} // namespace viz
} // namespace bsp
