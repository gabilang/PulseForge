// PulseForge — command-line interface.
//
// Usage:
//   pulseforge <command> [options]
//
// Commands:
//   demo                          Run a full demo on synthetic EEG/ECG/EMG/EOG.
//   eeg|ecg|emg|eog [options]     Analyse a channel (synthetic if no --in).
//   spectrum [options]            Print spectral summary of a signal.
//   filter   [options]            Filter a signal and write the result.
//   generate [options]           Generate a synthetic signal to --out.
//
// Common options:
//   --in <file.csv>   input CSV file (one column of samples)
//   --col <n>         column index (0-based, default 0)
//   --fs <hz>         sample rate (default depends on modality)
//   --dur <sec>       duration for synthetic signals
//   --out <file.csv>  output file
//
// filter options:  --ftype lowpass|highpass|bandpass|notch  --low <hz> --high <hz>
//                   --f0 <hz> --q <Q>
// generate options: --gtype ecg|eeg|emg|eog|sine  --freq <hz> --hr <bpm>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "bsp/Filters.hpp"
#include "bsp/Generators.hpp"
#include "bsp/Signal.hpp"
#include "bsp/Spectrum.hpp"
#include "bsp/Statistics.hpp"
#include "bsp/bio/ECG.hpp"
#include "bsp/bio/EEG.hpp"
#include "bsp/bio/EMG.hpp"
#include "bsp/bio/EOG.hpp"
#include "bsp/io/CsvIO.hpp"
#include "bsp/viz/AsciiPlot.hpp"
#include "bsp/viz/SvgPlot.hpp"

using namespace bsp;

namespace {

struct Args {
    std::map<std::string, std::string> opts;

    bool has(const std::string& k) const { return opts.count(k) > 0; }
    std::string str(const std::string& k, const std::string& def = "") const {
        auto it = opts.find(k);
        return it == opts.end() ? def : it->second;
    }
    double num(const std::string& k, double def) const {
        auto it = opts.find(k);
        return it == opts.end() ? def : std::stod(it->second);
    }
};

Args parseArgs(int argc, char** argv, int start) {
    Args a;
    for (int i = start; i < argc; ++i) {
        std::string tok = argv[i];
        if (tok.rfind("--", 0) == 0) {
            std::string key = tok.substr(2);
            if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0) {
                a.opts[key] = argv[++i];
            } else {
                a.opts[key] = "1";  // flag
            }
        }
    }
    return a;
}

void printHeader(const std::string& title) {
    std::cout << "\n========== " << title << " ==========\n";
}

// Resolve the signal to analyse: load from --in, else synthesise for `modality`.
Signal resolveSignal(const Args& a, const std::string& modality, double defaultFs) {
    double fs = a.num("fs", defaultFs);
    if (a.has("in")) {
        std::size_t col = static_cast<std::size_t>(a.num("col", 0));
        Signal s = io::loadSignal(a.str("in"), fs, col, modality);
        std::cout << "Loaded " << s.size() << " samples from " << a.str("in")
                  << " @ " << fs << " Hz\n";
        return s;
    }
    double dur = a.num("dur", 10.0);
    std::cout << "No --in given; synthesising " << dur << " s of " << modality
              << " @ " << fs << " Hz\n";
    if (modality == "ECG") return gen::syntheticECG(fs, dur, a.num("hr", 72.0));
    if (modality == "EEG") return gen::syntheticEEG(fs, dur);
    if (modality == "EMG") return gen::syntheticEMG(fs, dur);
    if (modality == "EOG") return gen::syntheticEOG(fs, dur, a.num("rate", 15.0));
    return gen::syntheticEEG(fs, dur);
}

void printBasicStats(const Signal& s) {
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  duration   : " << s.duration() << " s (" << s.size() << " samples)\n";
    std::cout << "  mean / RMS : " << stats::mean(s.samples) << " / " << stats::rms(s.samples) << "\n";
    std::cout << "  min / max  : " << stats::minValue(s.samples) << " / "
              << stats::maxValue(s.samples) << "\n";
    std::cout << "  std dev    : " << stats::stddev(s.samples) << "\n";
}

void printBanner();  // defined below; used by the demo and usage screens

// First `secs` seconds of a signal (for readable time-domain previews).
Signal headSeconds(const Signal& s, double secs) {
    std::size_t n = std::min(s.size(), static_cast<std::size_t>(secs * s.fs));
    return Signal(std::vector<double>(s.samples.begin(), s.samples.begin() + n), s.fs, s.label);
}

// Print an ASCII plot to the terminal when --ascii is given.
void showAscii(const Args& a, const std::string& text) {
    if (a.has("ascii")) std::cout << "\n" << text;
}

// Write an SVG to --svg when given.
void exportSvg(const Args& a, const std::string& svg) {
    if (!a.has("svg")) return;
    viz::writeSvg(a.str("svg"), svg);
    std::cout << "  wrote SVG -> " << a.str("svg") << "\n";
}

int cmdEEG(const Args& a) {
    Signal s = resolveSignal(a, "EEG", 256.0);
    printBasicStats(s);
    bio::EEGResult r = bio::analyzeEEG(s);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Band power (absolute | relative %):\n";
    auto row = [](const char* n, double abs, double rel) {
        std::cout << "    " << std::left << std::setw(7) << n << std::right
                  << std::setw(12) << abs << "   " << std::setw(6) << rel * 100.0 << " %\n";
    };
    row("delta", r.absolute.delta, r.relative.delta);
    row("theta", r.absolute.theta, r.relative.theta);
    row("alpha", r.absolute.alpha, r.relative.alpha);
    row("beta",  r.absolute.beta,  r.relative.beta);
    row("gamma", r.absolute.gamma, r.relative.gamma);
    std::cout << "  dominant freq      : " << r.dominantFreq << " Hz\n";
    std::cout << "  mean freq          : " << r.meanFreq << " Hz\n";
    std::cout << "  spectral edge 95%  : " << r.spectralEdge95 << " Hz\n";

    double fMax = a.num("fmax", 60.0);
    showAscii(a, viz::asciiSpectrum(r.psd, fMax, 78, 16, "EEG power spectrum"));
    if (a.has("svg")) {
        viz::SvgOptions o;
        o.title = "EEG power spectrum (Welch)";
        exportSvg(a, viz::svgSpectrum(r.psd, fMax, o));
    }
    return 0;
}

int cmdECG(const Args& a) {
    Signal s = resolveSignal(a, "ECG", 360.0);
    printBasicStats(s);
    bio::ECGResult r = bio::analyzeECG(s);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  R-peaks detected   : " << r.rPeaks.size() << "\n";
    std::cout << "  mean HR            : " << r.meanHR << " bpm\n";
    std::cout << "  HR range           : " << r.minHR << " - " << r.maxHR << " bpm\n";
    std::cout << "  SDNN  (HRV)        : " << r.sdnn << " ms\n";
    std::cout << "  RMSSD (HRV)        : " << r.rmssd << " ms\n";
    std::cout << "  pNN50              : " << r.pnn50 << " %\n";

    showAscii(a, viz::asciiWaveform(headSeconds(s, 6.0), 78, 16, "ECG (first 6 s)"));
    if (a.has("svg")) {
        viz::SvgOptions o;
        o.title = "ECG with detected R-peaks";
        exportSvg(a, viz::svgWaveform(s, o, r.rPeaks));
    }
    return 0;
}

int cmdEMG(const Args& a) {
    Signal s = resolveSignal(a, "EMG", 1000.0);
    printBasicStats(s);
    bio::EMGResult r = bio::analyzeEMG(s);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  RMS amplitude      : " << r.rms << "\n";
    std::cout << "  max envelope       : " << r.maxEnvelope << "\n";
    std::cout << "  activations        : " << r.activations.size() << "\n";
    for (std::size_t i = 0; i < r.activations.size(); ++i)
        std::cout << "      #" << (i + 1) << "  " << r.activations[i].startSec << " - "
                  << r.activations[i].endSec << " s  (" << r.activations[i].duration() << " s)\n";
    std::cout << std::setprecision(2);
    std::cout << "  mean frequency     : " << r.meanFrequency << " Hz\n";
    std::cout << "  median frequency   : " << r.medianFrequency << " Hz\n";

    Signal env(r.envelope, s.fs, "EMG envelope");
    showAscii(a, viz::asciiWaveform(env, 78, 14, "EMG RMS envelope"));
    if (a.has("svg")) {
        std::vector<std::size_t> onsets;
        for (const auto& act : r.activations) onsets.push_back(env.timeToIndex(act.startSec));
        viz::SvgOptions o;
        o.title = "EMG RMS envelope with activation onsets";
        o.lineColor = "#16a34a";
        exportSvg(a, viz::svgWaveform(env, o, onsets));
    }
    return 0;
}

int cmdEOG(const Args& a) {
    Signal s = resolveSignal(a, "EOG", 256.0);
    printBasicStats(s);
    bio::EOGResult r = bio::analyzeEOG(s);
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  blinks detected    : " << r.blinkCount << "\n";
    std::cout << "  blink rate         : " << r.blinkRatePerMin << " /min\n";
    std::cout << "  mean blink amp     : " << r.meanBlinkAmplitude << "\n";

    showAscii(a, viz::asciiWaveform(headSeconds(s, 20.0), 78, 14, "EOG (first 20 s)"));
    if (a.has("svg")) {
        viz::SvgOptions o;
        o.title = "EOG with detected blinks";
        o.lineColor = "#7c3aed";
        exportSvg(a, viz::svgWaveform(s, o, r.blinks));
    }
    return 0;
}

int cmdSpectrum(const Args& a) {
    Signal s = resolveSignal(a, "EEG", 256.0);
    printBasicStats(s);
    std::size_t seg = std::min<std::size_t>(s.size(), static_cast<std::size_t>(s.fs * 2.0));
    PSD psd = welch(s.samples, s.fs, seg, 0.5, WindowType::Hann);
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  freq resolution    : " << psd.resolution << " Hz\n";
    std::cout << "  dominant frequency : " << dominantFrequency(psd) << " Hz\n";
    std::cout << "  mean frequency     : " << meanFrequency(psd) << " Hz\n";
    std::cout << "  median frequency   : " << medianFrequency(psd) << " Hz\n";
    std::cout << "  spectral edge 95%  : " << spectralEdgeFrequency(psd, 0.95) << " Hz\n";
    std::cout << "  total power        : " << totalPower(psd) << "\n";

    double fMax = a.num("fmax", s.fs / 2.0);
    showAscii(a, viz::asciiSpectrum(psd, fMax, 78, 16, "power spectrum"));
    if (a.has("svg")) exportSvg(a, viz::svgSpectrum(psd, fMax));
    return 0;
}

int cmdFilter(const Args& a) {
    Signal s = resolveSignal(a, "EEG", 256.0);
    std::string type = a.str("ftype", "bandpass");
    std::vector<double> y;
    if (type == "lowpass")       y = filtfilt(s.samples, {biquadLowpass(a.num("high", 40.0), s.fs)});
    else if (type == "highpass") y = filtfilt(s.samples, {biquadHighpass(a.num("low", 1.0), s.fs)});
    else if (type == "notch")    y = notchFilter(s.samples, s.fs, a.num("f0", 50.0), a.num("q", 30.0));
    else                         y = bandpassZeroPhase(s.samples, s.fs, a.num("low", 1.0),
                                                       a.num("high", 40.0));
    std::cout << "Applied " << type << " filter.\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  input  RMS : " << stats::rms(s.samples) << "\n";
    std::cout << "  output RMS : " << stats::rms(y) << "\n";
    if (a.has("out")) {
        io::saveSignal(a.str("out"), Signal(y, s.fs, "filtered"));
        std::cout << "  wrote " << a.str("out") << "\n";
    }

    Signal filtered(y, s.fs, "filtered");
    if (a.has("ascii")) {
        std::cout << "\n" << viz::asciiWaveform(headSeconds(s, 4.0), 78, 12, "input (first 4 s)");
        std::cout << "\n" << viz::asciiWaveform(headSeconds(filtered, 4.0), 78, 12,
                                                "filtered (first 4 s)");
    }
    if (a.has("svg")) {
        viz::SvgOptions o;
        o.title = type + " filter output";
        exportSvg(a, viz::svgWaveform(filtered, o));
    }
    return 0;
}

int cmdPlot(const Args& a) {
    std::string modality = a.str("type", "EEG");
    std::transform(modality.begin(), modality.end(), modality.begin(), ::toupper);
    Signal s = resolveSignal(a, modality, a.num("fs", 256.0));
    printBasicStats(s);

    int w = static_cast<int>(a.num("width", 78));
    int h = static_cast<int>(a.num("height", 16));

    // Time-domain preview (windowed if long) and spectrum, both as ASCII.
    double previewSecs = a.num("secs", std::min(s.duration(), 10.0));
    std::cout << "\n" << viz::asciiWaveform(headSeconds(s, previewSecs), w, h,
                                            s.label + " (first " + std::to_string(
                                                static_cast<int>(previewSecs)) + " s)");
    std::size_t seg = std::min<std::size_t>(s.size(), static_cast<std::size_t>(s.fs * 2.0));
    PSD psd = welch(s.samples, s.fs, seg, 0.5, WindowType::Hann);
    double fMax = a.num("fmax", std::min(60.0, s.fs / 2.0));
    std::cout << "\n" << viz::asciiSpectrum(psd, fMax, w, h, "power spectrum");

    if (a.has("svg")) {
        viz::SvgOptions o;
        o.title = s.label + " waveform";
        exportSvg(a, viz::svgWaveform(s, o));
    }
    if (a.has("svg-spectrum")) {
        viz::writeSvg(a.str("svg-spectrum"), viz::svgSpectrum(psd, fMax));
        std::cout << "  wrote SVG -> " << a.str("svg-spectrum") << "\n";
    }
    return 0;
}

int cmdGenerate(const Args& a) {
    std::string type = a.str("gtype", "ecg");
    double fs = a.num("fs", 360.0);
    double dur = a.num("dur", 10.0);
    Signal s;
    if (type == "ecg")       s = gen::syntheticECG(fs, dur, a.num("hr", 72.0));
    else if (type == "eeg")  s = gen::syntheticEEG(fs, dur);
    else if (type == "emg")  s = gen::syntheticEMG(fs, dur);
    else if (type == "eog")  s = gen::syntheticEOG(fs, dur, a.num("rate", 15.0));
    else if (type == "sine") s = Signal(gen::sine(a.num("freq", 10.0), 1.0, fs, dur), fs, "sine");
    else { std::cerr << "Unknown --gtype " << type << "\n"; return 1; }

    std::string out = a.str("out", type + ".csv");
    io::saveSignal(out, s);
    std::cout << "Generated " << s.size() << " samples of " << type << " -> " << out << "\n";
    return 0;
}

int cmdDemo() {
    printBanner();
    std::cout << "Full demonstration on synthetic signals\n";

    printHeader("EEG (256 Hz, alpha-dominant)");
    {
        Signal s = gen::syntheticEEG(256.0, 12.0);
        bio::EEGResult r = bio::analyzeEEG(s);
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "  relative band power %:  delta " << r.relative.delta * 100
                  << "  theta " << r.relative.theta * 100
                  << "  alpha " << r.relative.alpha * 100
                  << "  beta " << r.relative.beta * 100
                  << "  gamma " << r.relative.gamma * 100 << "\n";
        std::cout << "  dominant frequency: " << std::setprecision(2) << r.dominantFreq
                  << " Hz (expected ~10)\n";
    }

    printHeader("ECG (360 Hz, 72 bpm)");
    {
        Signal s = gen::syntheticECG(360.0, 15.0, 72.0);
        bio::ECGResult r = bio::analyzeECG(s);
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "  R-peaks: " << r.rPeaks.size() << "   mean HR: " << r.meanHR
                  << " bpm (expected ~72)\n";
        std::cout << "  SDNN: " << r.sdnn << " ms   RMSSD: " << r.rmssd << " ms\n";
    }

    printHeader("EMG (1000 Hz, 3 bursts)");
    {
        Signal s = gen::syntheticEMG(1000.0, 10.0);
        bio::EMGResult r = bio::analyzeEMG(s);
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  activations detected: " << r.activations.size() << " (expected 3)\n";
        std::cout << "  median frequency: " << r.medianFrequency << " Hz\n";
    }

    printHeader("EOG (256 Hz, 15 blinks/min)");
    {
        Signal s = gen::syntheticEOG(256.0, 60.0, 15.0);
        bio::EOGResult r = bio::analyzeEOG(s);
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "  blinks: " << r.blinkCount << "   rate: " << r.blinkRatePerMin
                  << " /min (expected ~15)\n";
    }

    std::cout << "\nDemo complete.\n";
    return 0;
}

// Terminal logo: a stylised ECG/EEG pulse forged into the wordmark.
void printBanner() {
    std::cout <<
        "\n"
        "  ============================================================\n"
        "         /\\                  /\\                P U L S E F O R G E\n"
        "    ____/  \\____/\\_________/  \\______          \n"
        "              \\/                               Biosignal Processing Suite\n"
        "  ============================================================\n"
        "   EEG | ECG | EMG | EOG   -   signal, forged into insight\n"
        "\n";
}

void usage() {
    printBanner();
    std::cout <<
        "Usage: pulseforge <command> [options]\n\n"
        "Commands:\n"
        "  demo                 Run a full demo on synthetic EEG/ECG/EMG/EOG\n"
        "  eeg|ecg|emg|eog      Analyse a channel (synthetic if no --in)\n"
        "  spectrum             Spectral summary of a signal\n"
        "  filter               Filter a signal (--ftype lowpass|highpass|bandpass|notch)\n"
        "  generate             Write a synthetic signal (--gtype ecg|eeg|emg|eog|sine)\n"
        "  plot                 Visualise a signal (--type eeg|ecg|emg|eog)\n\n"
        "Common options: --in <csv> --col <n> --fs <hz> --dur <s> --out <csv>\n"
        "Visualisation:  --ascii          print an ASCII plot to the terminal\n"
        "                --svg <file>     export an SVG plot\n"
        "                --fmax <hz>      max frequency for spectra\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 0; }
    std::string cmd = argv[1];
    Args a = parseArgs(argc, argv, 2);

    try {
        if (cmd == "demo")     return cmdDemo();
        if (cmd == "eeg")      return cmdEEG(a);
        if (cmd == "ecg")      return cmdECG(a);
        if (cmd == "emg")      return cmdEMG(a);
        if (cmd == "eog")      return cmdEOG(a);
        if (cmd == "spectrum") return cmdSpectrum(a);
        if (cmd == "filter")   return cmdFilter(a);
        if (cmd == "generate") return cmdGenerate(a);
        if (cmd == "plot")     return cmdPlot(a);
        if (cmd == "help" || cmd == "--help" || cmd == "-h") { usage(); return 0; }
        std::cerr << "Unknown command: " << cmd << "\n\n";
        usage();
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}
