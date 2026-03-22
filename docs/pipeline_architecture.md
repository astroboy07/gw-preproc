# Pipeline Architecture

> Back to [docs README](README.md) |
> See also [preprocessing.md](preprocessing.md) for stage-by-stage details

## Goal

Produce a cleaned, narrowband complex baseband time series from raw
LIGO strain data. The output is suitable for any downstream analysis
that needs artifact-free, low-sample-rate data in a specific frequency
band.

## Data flow

```
  GWOSC HDF5 (4 kHz, 4096 s per file)
       │
       │  readligo.io.loaddata()
       │  -> strain (float64, 4096 Hz)
       │  -> channel_dict (DQ bitmasks at 1 Hz)
       │  -> gps_start (float)
       ▼
  ┌────────────────────────────────────────────────────┐
  │  Stage 1: Data-quality masking                     │
  │  apply_dq_mask(strain, channel_dict, fs)           │
  │                                                    │
  │  - Combine DATA and CBC_CAT1 flags into a 1 Hz     │
  │    boolean mask                                    │
  │  - Stretch to 4096 Hz, zero bad samples            │
  │                                                    │
  │  Output: clean_strain (4096 Hz) + dq_mask_1hz      │
  └────────────────────────┬───────────────────────────┘
                           │
                           ▼
  ┌────────────────────────────────────────────────────┐
  │  Stage 2: Heterodyne and decimation                │
  │  heterodyne_downsample(clean_strain, fs, f_band)   │
  │                                                    │
  │  - Multiply by exp(-2*pi*i*f_h*t) to shift band    │
  │    centre to DC                                    │
  │  - Multi-stage FIR decimation (each factor <= 13)  │
  │    until sample rate = bandwidth                   │
  │                                                    │
  │  Output: narrowband (complex, fs_new Hz)           │
  └────────────────────────┬───────────────────────────┘
                           │
                           ▼
  ┌────────────────────────────────────────────────────┐
  │  Stage 3: STFT-based transient cleaning            │
  │  clean_narrowband(narrowband, dq_mask_1hz, fs_new) │
  │                                                    │
  │  - Two-sided STFT (Hann, 50% overlap)              │
  │  - Normalise power by 39.4th percentile of         │
  │    DQ-good columns (-> chi-squared(2) under H0)   │
  │  - Flag columns with > bound bins above            │
  │    chi-squared threshold at FAP = fap              │
  │  - Zero flagged columns, reconstruct via ISTFT     │
  │                                                    │
  │  Output: cleaned (complex) + sample_mask (bool)    │
  └────────────────────────────────────────────────────┘
```

## Entry point

All three stages are chained by `process_band()`:

```python
gps_times, cleaned, sample_mask, fs_new = process_band(
    strain, channel_dict, gps_start, f_band=(100, 101)
)
```

## Module responsibilities

| Module | Responsibility |
|---|---|
| `readligo.io` | Read GWOSC HDF5 files, extract strain + DQ bitmasks + metadata |
| `readligo.preprocessing` | Per-file narrowband cleaning pipeline (all three stages above) |

## Performance notes

- **Multi-stage decimation** keeps each `scipy.signal.decimate` call
  at factor <= 13 (FIR mode) for numerical stability with complex data.
- **STFT cleaning** operates at the decimated sample rate (e.g. 1 Hz
  for a 1 Hz band), so it is very fast relative to the raw data size.
- The **39.4th percentile normalisation** is chosen because for an
  exponentially distributed variable (|Z|^2 under Gaussian noise),
  the 39.4th percentile equals the mean, mapping the normalised power
  to chi-squared(2).
