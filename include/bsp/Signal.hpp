#pragma once
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace bsp {

// A uniformly-sampled, single-channel signal: raw samples plus a sample rate.
class Signal {
public:
    std::vector<double> samples;
    double fs = 1.0;           // sample rate in Hz
    std::string label;         // human readable channel name

    Signal() = default;
    Signal(std::vector<double> s, double rate, std::string lbl = "")
        : samples(std::move(s)), fs(rate), label(std::move(lbl)) {}

    std::size_t size() const { return samples.size(); }
    bool empty() const { return samples.empty(); }
    double duration() const { return fs > 0.0 ? static_cast<double>(samples.size()) / fs : 0.0; }

    double& operator[](std::size_t i) { return samples[i]; }
    double operator[](std::size_t i) const { return samples[i]; }

    // Index of the sample nearest to time t (seconds).
    std::size_t timeToIndex(double t) const {
        long idx = static_cast<long>(t * fs + 0.5);
        if (idx < 0) idx = 0;
        if (static_cast<std::size_t>(idx) >= samples.size() && !samples.empty())
            idx = static_cast<long>(samples.size()) - 1;
        return static_cast<std::size_t>(idx);
    }
};

} // namespace bsp
