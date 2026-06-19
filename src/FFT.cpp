#include "bsp/FFT.hpp"

#include <cmath>

namespace bsp {

static constexpr double kPi = 3.14159265358979323846;

std::size_t nextPowerOfTwo(std::size_t n) {
    std::size_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

void fft(CArray& a) {
    const std::size_t n = a.size();
    if (n < 2) return;

    // Bit-reversal permutation.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // Danielson-Lanczos butterflies.
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const double ang = -2.0 * kPi / static_cast<double>(len);
        const Complex wlen(std::cos(ang), std::sin(ang));
        for (std::size_t i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            for (std::size_t k = 0; k < len / 2; ++k) {
                Complex u = a[i + k];
                Complex v = a[i + k + len / 2] * w;
                a[i + k] = u + v;
                a[i + k + len / 2] = u - v;
                w *= wlen;
            }
        }
    }
}

void ifft(CArray& a) {
    for (auto& c : a) c = std::conj(c);
    fft(a);
    const double inv = 1.0 / static_cast<double>(a.size());
    for (auto& c : a) c = std::conj(c) * inv;
}

CArray rfft(const std::vector<double>& x) {
    const std::size_t n = nextPowerOfTwo(x.size());
    CArray a(n, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < x.size(); ++i) a[i] = Complex(x[i], 0.0);
    fft(a);
    return a;
}

} // namespace bsp
