// Lightweight self-contained test runner (no external framework).
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "bsp/Filters.hpp"
#include "bsp/Generators.hpp"
#include "bsp/Spectrum.hpp"
#include "bsp/Statistics.hpp"
#include "bsp/bio/ECG.hpp"
#include "bsp/bio/EEG.hpp"
#include "bsp/bio/EMG.hpp"
#include "bsp/bio/EOG.hpp"
#include "bsp/viz/AsciiPlot.hpp"
#include "bsp/viz/SvgPlot.hpp"

using namespace bsp;

static int g_failures = 0;
static int g_checks = 0;

static void check(bool cond, const std::string& msg) {
    ++g_checks;
    if (!cond) {
        ++g_failures;
        std::cerr << "  [FAIL] " << msg << "\n";
    } else {
        std::cout << "  [ok]   " << msg << "\n";
    }
}

static void checkNear(double got, double want, double tol, const std::string& msg) {
    bool ok = std::fabs(got - want) <= tol;
    check(ok, msg + " (got " + std::to_string(got) + ", want " + std::to_string(want) + ")");
}

static void testStatistics() {
    std::cout << "Statistics\n";
    std::vector<double> x = {1, 2, 3, 4, 5};
    checkNear(stats::mean(x), 3.0, 1e-9, "mean");
    checkNear(stats::median(x), 3.0, 1e-9, "median");
    checkNear(stats::stddev(x), std::sqrt(2.5), 1e-9, "stddev (sample)");
    checkNear(stats::rms(std::vector<double>{3, 4}), std::sqrt(12.5), 1e-9, "rms");
}

static void testSpectrum() {
    std::cout << "Spectrum / FFT\n";
    double fs = 256.0;
    auto x = gen::sine(10.0, 1.0, fs, 4.0);
    PSD psd = welch(x, fs, 512, 0.5, WindowType::Hann);
    checkNear(dominantFrequency(psd), 10.0, 1.0, "dominant freq of 10 Hz sine");

    // A two-tone signal: most power should sit below 30 Hz.
    auto a = gen::sine(5.0, 1.0, fs, 4.0);
    auto b = gen::sine(20.0, 1.0, fs, 4.0);
    std::vector<double> mix(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) mix[i] = a[i] + b[i];
    PSD pm = welch(mix, fs, 512, 0.5, WindowType::Hann);
    double low = bandPower(pm, 0, 30), high = bandPower(pm, 30, 128);
    check(low > high * 10.0, "two-tone power concentrated below 30 Hz");
}

static void testFilters() {
    std::cout << "Filters\n";
    double fs = 256.0;
    auto lowTone = gen::sine(5.0, 1.0, fs, 4.0);
    auto highTone = gen::sine(60.0, 1.0, fs, 4.0);
    std::vector<double> mix(lowTone.size());
    for (std::size_t i = 0; i < mix.size(); ++i) mix[i] = lowTone[i] + highTone[i];

    auto lp = filtfilt(mix, {biquadLowpass(15.0, fs)});
    PSD before = welch(mix, fs, 512, 0.5, WindowType::Hann);
    PSD after = welch(lp, fs, 512, 0.5, WindowType::Hann);
    double hiBefore = bandPower(before, 40, 128);
    double hiAfter = bandPower(after, 40, 128);
    check(hiAfter < hiBefore * 0.05, "lowpass attenuates 60 Hz tone by >20x");

    // Notch should remove a 50 Hz interferer.
    auto sig = gen::sine(50.0, 1.0, fs, 4.0);
    auto notched = notchFilter(sig, fs, 50.0, 30.0);
    check(stats::rms(notched) < stats::rms(sig) * 0.3, "50 Hz notch attenuates interferer");

    // FIR lowpass should have ~unity DC gain.
    auto h = designLowpassFIR(101, 20.0, fs);
    double dc = 0.0;
    for (double v : h) dc += v;
    checkNear(dc, 1.0, 1e-6, "FIR lowpass unity DC gain");
}

static void testECG() {
    std::cout << "ECG\n";
    double fs = 360.0, dur = 20.0, hr = 72.0;
    Signal s = gen::syntheticECG(fs, dur, hr);
    bio::ECGResult r = bio::analyzeECG(s);
    int expected = static_cast<int>(dur * hr / 60.0);
    check(std::abs(static_cast<int>(r.rPeaks.size()) - expected) <= 2,
          "R-peak count near expected (" + std::to_string(expected) + ")");
    checkNear(r.meanHR, hr, 5.0, "mean heart rate");
}

static void testEEG() {
    std::cout << "EEG\n";
    Signal s = gen::syntheticEEG(256.0, 12.0);
    bio::EEGResult r = bio::analyzeEEG(s);
    check(r.relative.alpha > r.relative.gamma, "alpha power exceeds gamma (alpha-dominant)");
    checkNear(r.dominantFreq, 10.0, 1.5, "dominant frequency near alpha (10 Hz)");
}

static void testEMG() {
    std::cout << "EMG\n";
    Signal s = gen::syntheticEMG(1000.0, 10.0);
    bio::EMGResult r = bio::analyzeEMG(s);
    check(r.activations.size() == 3, "detects 3 muscle activations");
}

static void testEOG() {
    std::cout << "EOG\n";
    Signal s = gen::syntheticEOG(256.0, 60.0, 15.0);
    bio::EOGResult r = bio::analyzeEOG(s);
    check(std::abs(r.blinkCount - 15) <= 2, "blink count near 15");
}

static void testViz() {
    std::cout << "Visualization\n";
    Signal s = gen::syntheticECG(360.0, 4.0, 72.0);
    std::string ascii = viz::asciiWaveform(s, 60, 12, "ecg");
    check(ascii.find('\n') != std::string::npos && ascii.size() > 100, "ascii waveform non-empty");

    bio::ECGResult r = bio::analyzeECG(s);
    std::string svg = viz::svgWaveform(s, {}, r.rPeaks);
    check(svg.rfind("<svg", 0) == 0, "svg starts with <svg");
    check(svg.find("</svg>") != std::string::npos, "svg closed properly");
    // one <circle per detected R-peak marker.
    std::size_t circles = 0, pos = 0;
    while ((pos = svg.find("<circle", pos)) != std::string::npos) { ++circles; pos += 7; }
    check(circles == r.rPeaks.size(), "svg has one marker per R-peak");

    PSD psd = welch(s.samples, s.fs, 512, 0.5, WindowType::Hann);
    std::string spec = viz::svgSpectrum(psd, 60.0);
    check(spec.find("<path") != std::string::npos, "svg spectrum has a path");
}

int main() {
    std::cout << "Running BioSignalProcessor tests\n\n";
    testStatistics();
    testSpectrum();
    testFilters();
    testECG();
    testEEG();
    testEMG();
    testEOG();
    testViz();

    std::cout << "\n" << (g_checks - g_failures) << "/" << g_checks << " checks passed.\n";
    if (g_failures > 0) {
        std::cout << g_failures << " FAILURES\n";
        return 1;
    }
    std::cout << "ALL TESTS PASSED\n";
    return 0;
}
