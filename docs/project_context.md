# Project Context

Last updated: 2026-03-22

## What this project does

Reads and cleans **LIGO strain data** (O3a run, GWOSC 4 kHz HDF5 format)
to produce narrowband complex baseband time series suitable for downstream
gravitational-wave analyses. The cleaning pipeline removes data-quality
flagged segments and transient artifacts (glitches) that survive the DQ
flags.

## Key libraries

| Library | Role |
|---|---|
| `numpy`, `scipy` | Core numerics, FFTs, signal processing |
| `h5py` | Reading GWOSC HDF5 strain files |
| `matplotlib` | Diagnostic plots |

## Repository layout

```
readligo/
├── readligo/              # Installable Python package
│   ├── __init__.py        # Public API re-exports
│   ├── io.py              # GWOSC HDF5 reader (loaddata, read_hdf5)
│   └── preprocessing.py   # DQ masking, heterodyne, STFT cleaning
├── notebooks/
│   └── data_preprocessing_v3.ipynb   # Interactive exploration
├── docs/
│   ├── README.md                     # Documentation hub (start here)
│   ├── project_context.md            # You are here
│   ├── pipeline_architecture.md      # Data flow and module design
│   ├── preprocessing.md              # Preprocessing deep-dive
│   └── rendered/                     # PDFs and TeX sources
│       ├── cleaning_summary.pdf
│       └── cleaning_summary.tex
├── pyproject.toml
└── pyrightconfig.json
```

## Key files and their purposes

| File | Purpose |
|---|---|
| `readligo/io.py` | Reads GWOSC HDF5 files. Extracts strain, DQ bitmasks, injection masks, GPS metadata. Returns typed numpy arrays. |
| `readligo/preprocessing.py` | Full single-file preprocessing: DQ masking, heterodyne to baseband, multi-stage decimation, STFT-based glitch removal with chi-squared test. |
| `readligo/__init__.py` | Re-exports `loaddata`, `read_hdf5`, `process_band`, `apply_dq_mask`, `heterodyne_downsample`, `clean_narrowband`, `plot_cleaning`. |
| `notebooks/data_preprocessing_v3.ipynb` | Step-by-step walkthrough of the preprocessing pipeline on a single O3a file. |

## Recent achievements

- **Narrowband preprocessing pipeline**: end-to-end from raw 4 kHz
  strain to cleaned 1 Hz complex baseband, with per-sample validity
  masks and diagnostic plotting.
- **Package structure**: installable via `pip install -e .`, clean type
  annotations, Pylance-clean.

## Conventions

- Python 3.11, type-annotated, Pylance-clean
- `numpy` arrays everywhere — no pandas in the core pipeline
- GPS times as `float64` (seconds since 1980-01-06T00:00:00 UTC)
- Frequency bands specified as `(f_min, f_max)` tuples in Hz
- Complex baseband convention: `exp(-2 pi i f_h t)` heterodyne
