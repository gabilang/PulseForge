#pragma once
#include <vector>

#include "bsp/Windows.hpp"

namespace bsp {

// Second-order IIR section (transposed direct form II), RBJ "audio EQ
// cookbook" coefficients. a0 is normalised to 1.
struct Biquad {
    double b0 = 1.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0;
    double z1 = 0.0, z2 = 0.0;  // state

    double process(double x);
    void reset() { z1 = z2 = 0.0; }
};

Biquad biquadLowpass(double fc, double fs, double Q = 0.70710678);
Biquad biquadHighpass(double fc, double fs, double Q = 0.70710678);
Biquad biquadBandpass(double fc, double fs, double Q);
Biquad biquadNotch(double fc, double fs, double Q);

// Causal application of a single biquad (state runs forward over x).
std::vector<double> applyBiquad(Biquad bq, const std::vector<double>& x);

// Zero-phase forward-backward filtering through a cascade of biquads,
// with reflection padding to suppress edge transients.
std::vector<double> filtfilt(const std::vector<double>& x, std::vector<Biquad> chain);

// Windowed-sinc linear-phase FIR design. numTaps is forced odd.
std::vector<double> designLowpassFIR(int numTaps, double fc, double fs,
                                     WindowType w = WindowType::Hamming);
std::vector<double> designHighpassFIR(int numTaps, double fc, double fs,
                                      WindowType w = WindowType::Hamming);
std::vector<double> designBandpassFIR(int numTaps, double f1, double f2, double fs,
                                      WindowType w = WindowType::Hamming);

// Apply a symmetric FIR with group-delay compensation -> same length, zero phase.
std::vector<double> applyFIR(const std::vector<double>& h, const std::vector<double>& x);

// Simple sliding moving-average smoother (odd window recommended).
std::vector<double> movingAverage(const std::vector<double>& x, int window);

// Convenience powerline notch (e.g. 50 or 60 Hz) via zero-phase biquad.
std::vector<double> notchFilter(const std::vector<double>& x, double fs, double f0, double Q = 30.0);

// Convenience zero-phase bandpass via highpass+lowpass biquad cascade.
std::vector<double> bandpassZeroPhase(const std::vector<double>& x, double fs,
                                      double low, double high, double Q = 0.70710678);

} // namespace bsp
