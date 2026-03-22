# readligo

Python package for reading and preprocessing LIGO gravitational-wave strain data (GWOSC HDF5 format).

## Installation

```bash
pip install -e .
```

## Quick start

```python
from readligo import loaddata, process_band

# Load strain and data-quality channels from an HDF5 file
strain, time, channel_dict = loaddata("H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5")
gps_start = time["start"]

# Preprocess a 1 Hz band around 100 Hz
gps_times, cleaned, sample_mask, fs_new = process_band(
    strain, channel_dict, gps_start, f_band=(100, 101)
)
```

## Documentation

See [`docs/README.md`](docs/README.md) for the full documentation hub, including:

- [Project context](docs/project_context.md) — status, architecture, conventions
- [Pipeline architecture](docs/pipeline_architecture.md) — algorithm overview
- [Preprocessing](docs/preprocessing.md) — detailed cleaning pipeline

## Package layout

```
readligo/
├── readligo/
│   ├── __init__.py
│   ├── io.py              # HDF5 reader (loaddata, read_hdf5)
│   └── preprocessing.py   # Signal processing pipeline
├── notebooks/
│   └── data_preprocessing_v3.ipynb
├── docs/
│   ├── README.md           # Documentation hub
│   ├── project_context.md
│   ├── pipeline_architecture.md
│   ├── preprocessing.md
│   └── rendered/           # PDFs and TeX sources
├── pyproject.toml
└── pyrightconfig.json
```
