#pragma once
#include <cstdint>
#include <vector>

#include "bsp/Signal.hpp"

namespace bsp {
namespace gen {

// Basic waveforms.
std::vector<double> sine(double freqHz, double amp, double fs, double durationSec, double phase = 0.0);
std::vector<double> whiteNoise(double sigma, double fs, double durationSec, std::uint32_t seed = 1234);

// Synthetic biosignals so the app is fully runnable without external data.
// hr is heart rate in beats per minute.
Signal syntheticECG(double fs, double durationSec, double hr = 60.0, double noise = 0.01,
                    std::uint32_t seed = 7);
// Sum of band oscillations (delta..gamma) plus broadband noise.
Signal syntheticEEG(double fs, double durationSec, double noise = 5.0, std::uint32_t seed = 11);
// Quiet baseline punctuated by burst activations (one per second window listed).
Signal syntheticEMG(double fs, double durationSec, std::uint32_t seed = 23);
// Slow drift plus discrete blink bumps; blinkRate is blinks per minute.
Signal syntheticEOG(double fs, double durationSec, double blinkRate = 15.0, std::uint32_t seed = 41);

} // namespace gen
} // namespace bsp
