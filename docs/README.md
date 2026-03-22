# readligo — Documentation Hub

Package for reading and cleaning LIGO/GWOSC strain data. Reads HDF5
files, applies data-quality masks, heterodynes to a narrowband baseband,
and removes transient artifacts via STFT-based cleaning.

## Quick navigation

| Document | Purpose |
|---|---|
| [project_context.md](project_context.md) | Current status, key files, conventions |
| [pipeline_architecture.md](pipeline_architecture.md) | Preprocessing data flow and module responsibilities |
| [preprocessing.md](preprocessing.md) | Detailed preprocessing stages (DQ, heterodyne, STFT cleaning) |
| [rendered/](rendered/) | PDFs and TeX source (cleaning_summary) |

## Pipeline at a glance

```
  GWOSC HDF5 files (4 kHz strain, 4096-second segments)
        │
        ▼
  ┌─────────────┐
  │  readligo.io │   loaddata() / read_hdf5()
  └──────┬──────┘
         │  strain, channel_dict, gps_start
         ▼
  ┌──────────────────────┐
  │  preprocessing       │
  │  ├ apply_dq_mask     │   zero bad-DQ samples
  │  ├ heterodyne_       │   shift band to baseband,
  │  │  downsample       │   decimate to 1 Hz
  │  └ clean_narrowband  │   STFT glitch removal
  └──────┬───────────────┘
         │
         ▼
  cleaned complex baseband + sample_mask
```

## Recent progress

- Package structure set up (`readligo` installable via `pip install -e .`)
- I/O module reads GWOSC HDF5 files with full type safety
- Narrowband preprocessing pipeline operational (DQ masking, heterodyne
  decimation, STFT-based transient cleaning with chi-squared test)

## For new sessions

If you are an AI assistant starting a new session on this project:

1. **Read this file first** for high-level orientation.
2. **Read [project_context.md](project_context.md)** for current status
   and what to prioritize.
3. **Read [pipeline_architecture.md](pipeline_architecture.md)** for the
   data flow and how modules connect.
4. **Read [preprocessing.md](preprocessing.md)** for detailed stage-by-stage
   documentation.
5. Source code lives in `readligo/` — `io.py` (I/O) and
   `preprocessing.py` (signal processing).
6. Exploratory work is in `notebooks/`.
