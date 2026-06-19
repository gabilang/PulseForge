#pragma once
#include <cstddef>
#include <vector>

#include "bsp/Signal.hpp"

namespace bsp {
namespace bio {

struct EMGActivation {
    double startSec = 0;
    double endSec = 0;
    double duration() const { return endSec - startSec; }
};

struct EMGResult {
    std::vector<double> envelope;             // RMS envelope (same length as input)
    std::vector<EMGActivation> activations;   // detected muscle-activation bursts
    double rms = 0;                           // overall RMS amplitude
    double maxEnvelope = 0;
    double meanFrequency = 0;                 // spectral mean frequency (Hz)
    double medianFrequency = 0;               // spectral median frequency (Hz, fatigue index)
};

EMGResult analyzeEMG(const Signal& sig);

} // namespace bio
} // namespace bsp
