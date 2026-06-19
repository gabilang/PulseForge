#pragma once
#include "bsp/Signal.hpp"
#include "bsp/Spectrum.hpp"

namespace bsp {
namespace bio {

// Classic EEG frequency bands (Hz).
struct EEGBands {
    double deltaLo = 0.5, deltaHi = 4.0;
    double thetaLo = 4.0, thetaHi = 8.0;
    double alphaLo = 8.0, alphaHi = 13.0;
    double betaLo = 13.0, betaHi = 30.0;
    double gammaLo = 30.0, gammaHi = 100.0;
};

struct BandPower {
    double delta = 0, theta = 0, alpha = 0, beta = 0, gamma = 0, total = 0;
};

struct EEGResult {
    BandPower absolute;        // power in each band
    BandPower relative;        // fraction of total power in each band
    double dominantFreq = 0;   // peak frequency (Hz)
    double spectralEdge95 = 0; // 95% spectral edge frequency (Hz)
    double meanFreq = 0;
    PSD psd;                   // underlying Welch PSD
};

EEGResult analyzeEEG(const Signal& sig, EEGBands bands = {});

} // namespace bio
} // namespace bsp
