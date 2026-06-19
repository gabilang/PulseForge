#pragma once
#include <cstddef>
#include <vector>

#include "bsp/Signal.hpp"

namespace bsp {
namespace bio {

struct ECGResult {
    std::vector<std::size_t> rPeaks;     // R-peak sample indices
    std::vector<double> rrIntervals;     // successive R-R intervals (seconds)
    double meanHR = 0;                   // mean heart rate (bpm)
    double minHR = 0;
    double maxHR = 0;
    double sdnn = 0;                     // HRV: std of R-R intervals (ms)
    double rmssd = 0;                    // HRV: RMS of successive differences (ms)
    double pnn50 = 0;                    // HRV: % of |ΔRR| > 50 ms
};

// Simplified Pan–Tompkins QRS detection plus heart-rate / HRV metrics.
ECGResult analyzeECG(const Signal& sig);

} // namespace bio
} // namespace bsp
