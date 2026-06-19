<p align="center">
  <img src="assets/pulseforge-logo.svg" alt="PulseForge" width="380"/>
</p>

# PulseForge — Biosignal Processing Suite

*Signal, forged into insight.*

A C++17 toolkit and command-line application for **digital signal processing of
biosignals** — EEG, ECG, EMG and EOG. It bundles a reusable DSP core (FFT,
filters, spectral estimation, statistics) with modality-specific analysis
pipelines, plus synthetic signal generators so the whole thing runs with zero
external data.

## Features

### DSP core (`bsp` library)
- **FFT** — iterative radix-2 Cooley–Tukey, forward/inverse, real-input helper.
- **Spectral estimation** — Welch's method and periodograms, band power
  (trapezoidal integration), dominant / mean / median / spectral-edge frequency.
- **Filters**
  - Zero-phase IIR (`filtfilt`) via RBJ biquads: low-pass, high-pass, band-pass, notch.
  - Linear-phase windowed-sinc FIR design (low/high/band-pass) with group-delay compensation.
  - Power-line notch (50/60 Hz) and moving-average helpers.
- **Windows** — Hann, Hamming, Blackman, rectangular.
- **Statistics** — mean, variance, RMS, median, skewness, excess kurtosis, zero-crossing rate.
- **Generators** — sine, Gaussian white noise, and realistic synthetic EEG/ECG/EMG/EOG.
- **Visualizers** — terminal ASCII plots (waveform & spectrum) and standalone SVG
  export, with biosignal-aware overlays (ECG R-peaks, EOG blinks, EMG activation onsets).

### Biosignal pipelines
| Modality | Analysis |
|----------|----------|
| **EEG** | Welch PSD, absolute & relative band power (delta/theta/alpha/beta/gamma), dominant & spectral-edge frequency |
| **ECG** | Simplified Pan–Tompkins QRS detection, heart rate, HRV (SDNN, RMSSD, pNN50) |
| **EMG** | Band-pass + RMS envelope, activation/onset detection, mean & median frequency (fatigue indices) |
| **EOG** | Blink detection, blink rate, mean blink amplitude |

## Build

Requires a C++17 compiler. Two equivalent build paths:

```bash
# Plain Makefile (no extra tooling)
make            # builds the app and runs the tests
make app        # just the CLI -> build/pulseforge
make test       # build + run the test suite
make run        # build + run the demo

# Or with CMake
cmake -S . -B build && cmake --build build
ctest --test-dir build
```

> On some macOS Command Line Tools installs the toolchain-local libc++ headers
> are incomplete; the Makefile automatically falls back to the SDK's headers
> (`xcrun --show-sdk-path`). With CMake, if you hit a `'complex' file not found`
> error, configure with
> `-DCMAKE_CXX_FLAGS="-isystem $(xcrun --show-sdk-path)/usr/include/c++/v1"`.

## Usage

```bash
pulseforge demo                         # full demo on synthetic signals

pulseforge eeg                          # analyse synthetic EEG (256 Hz)
pulseforge ecg --in ecg.csv --fs 360    # analyse a recording from CSV
pulseforge emg --fs 1000 --dur 10
pulseforge eog --rate 20

pulseforge spectrum --in eeg.csv --fs 256
pulseforge filter --ftype notch --f0 50 --in eeg.csv --fs 256 --out clean.csv
pulseforge generate --gtype ecg --fs 360 --dur 30 --hr 75 --out ecg.csv
pulseforge plot --type ecg --fs 360                 # ASCII waveform + spectrum
```

**Common options:** `--in <csv>` `--col <n>` `--fs <hz>` `--dur <sec>` `--out <csv>`.
If `--in` is omitted, the relevant modality is synthesised so every command is runnable.

### Visualization
Add `--ascii` to print plots in the terminal, or `--svg <file>` to export a
standalone SVG (opens in any browser). The plot is tailored to the command:

```bash
pulseforge ecg --ascii                       # ECG waveform in the terminal
pulseforge ecg --svg ecg.svg                  # SVG waveform + R-peak markers
pulseforge eeg --ascii --fmax 40              # EEG power spectrum (0-40 Hz)
pulseforge eog --svg eog.svg                  # EOG trace + blink markers
pulseforge emg --svg emg.svg                  # RMS envelope + activation onsets
pulseforge plot --type eeg --svg w.svg --svg-spectrum s.svg
```

Options: `--ascii`, `--svg <file>`, `--fmax <hz>` (spectrum limit), `--width`/`--height`
(ASCII size, `plot` command).

### Input format
A signal CSV is **one sample per line** (an optional header row is tolerated).
The sample rate is supplied with `--fs`; the time axis is implicit. Multi-column
files are supported via `--col <n>` (0-based).

## Project layout
```
include/bsp/        public headers (Signal, FFT, Filters, Spectrum, bio/*, viz/*)
src/                library implementation (incl. src/viz ASCII + SVG plotters)
app/main.cpp        command-line interface
tests/test_main.cpp self-contained test runner
CMakeLists.txt      CMake build
Makefile            plain Make build
```

## Testing
`make test` runs a lightweight self-contained suite that validates the FFT/PSD
(peak frequency), filters (band attenuation, unity DC gain, notch), and each
pulseforge pipeline against synthetic ground truth (e.g. ECG beat count at a known
heart rate, alpha-dominant EEG, EMG burst count, EOG blink count).
# PulseForge
