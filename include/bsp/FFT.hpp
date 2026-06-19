#pragma once
#include <complex>
#include <cstddef>
#include <vector>

namespace bsp {

using Complex = std::complex<double>;
using CArray = std::vector<Complex>;

// Smallest power of two >= n.
std::size_t nextPowerOfTwo(std::size_t n);

// In-place iterative radix-2 Cooley-Tukey FFT. Size must be a power of two.
void fft(CArray& a);

// In-place inverse FFT (normalised by 1/N).
void ifft(CArray& a);

// Real-input FFT. The input is zero-padded to the next power of two.
CArray rfft(const std::vector<double>& x);

} // namespace bsp
