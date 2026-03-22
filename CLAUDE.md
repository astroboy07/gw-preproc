# CLAUDE.md

Read `docs/README.md` first for project orientation, then
`docs/project_context.md` for current status and conventions.
Only go deeper into the specific docs below when the task requires it.

## Documentation map

| File | Read when... |
|---|---|
| `docs/README.md` | Always — start here |
| `docs/project_context.md` | Always — current state, libraries, conventions |
| `docs/pipeline_architecture.md` | Working on how the three stages connect |
| `docs/preprocessing.md` | Working on DQ masking, heterodyne, or STFT cleaning |

## Source code

| File | Purpose |
|---|---|
| `readligo/io.py` | GWOSC HDF5 reader — `loaddata()`, `read_hdf5()` |
| `readligo/preprocessing.py` | Cleaning pipeline — `process_band()` and its sub-functions |
| `readligo/__init__.py` | Public API re-exports |

## Notebooks

`notebooks/data_preprocessing.ipynb` — interactive walkthrough of the
full preprocessing pipeline on a single O3a file.

## Key conventions

- Python 3.11, fully type-annotated, Pylance-clean
- `numpy` arrays throughout; GPS times as `float64`
- Frequency bands as `(f_min, f_max)` tuples in Hz
- Complex baseband uses `exp(-2 pi i f_h t)` heterodyne convention
- `loaddata()` raises `ValueError` on empty or non-HDF5 files (never returns `None`)
