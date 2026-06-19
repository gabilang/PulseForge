#pragma once
#include <cstddef>
#include <vector>

#include "bsp/Signal.hpp"

namespace bsp {
namespace bio {

struct EOGResult {
    std::vector<std::size_t> blinks;   // blink peak sample indices
    int blinkCount = 0;
    double blinkRatePerMin = 0;
    double meanBlinkAmplitude = 0;
};

// Detect eye blinks (large transient deflections) in an EOG channel.
EOGResult analyzeEOG(const Signal& sig);

} // namespace bio
} // namespace bsp
