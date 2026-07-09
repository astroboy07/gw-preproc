import h5py
import numpy as np
from pathlib import Path
from rich.progress import track
from .io import loaddata
from .preprocessing import process_band


DATA_DIRNAME = "/Users/saifa/Library/CloudStorage/OneDrive-weizmann.ac.il/O3a_marvin/"
OUT_DIRNAME = "/Users/saifa/Library/CloudStorage/OneDrive-weizmann.ac.il/O3a_marvin_processed"

def save_hdf5_file(
    out_dir: str | Path,
    source_filename: str,
    f_band: tuple[float, float],
    gps_times: np.ndarray,
    cleaned: np.ndarray,
    sample_mask: np.ndarray,
    fs_new: float,
    breakdown: dict[str, int | float] | None,
) -> None:
    out_path = Path(out_dir) / (
        Path(source_filename).stem + f"_f_band_{int(f_band[0])}-{int(f_band[1])}Hz.hdf5"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    
    with h5py.File(out_path, "w") as f:
        # ── metadata ──────────────────────────────────────────────
        f.attrs["source_file"]  = source_filename
        f.attrs["f_band"]       = list(f_band)   # [f_min, f_max] Hz
        f.attrs["fs_new"]       = fs_new

        # ── breakdown stats ────────────────────────────────────────
        if breakdown is not None:
            stats = f.require_group("stats")
            for key, val in breakdown.items():
                stats.attrs[key] = val

        # ── time series ───────────────────────────────────────────
        f.create_dataset("gps_times",   data=gps_times,               compression="gzip")
        f.create_dataset("cleaned_re",  data=np.real(cleaned),            compression="gzip")
        f.create_dataset("cleaned_im",  data=np.imag(cleaned),            compression="gzip")
        f.create_dataset("sample_mask", data=sample_mask.astype(np.uint8), compression="gzip")


def process_ligo_files(
    f_band: tuple[float, float],
    nperseg: int = 64,
    fap: float = 0.001,
    bound: int = 10,
    data_dirname: str = DATA_DIRNAME,
    process_dirname: str = OUT_DIRNAME,
) -> None:
    hdf5_files = list(Path(data_dirname).glob("*.hdf5"))
    hdf5_files_sorted = sorted(
        hdf5_files,
        key=lambda p: int(p.stem.split('-')[-2])
    )
    for files in track(hdf5_files_sorted, description="Processing files"):
        strain, time, channel_dict = loaddata(data_dirname + files.name)
        assert isinstance(time, np.ndarray), "expected time vector"
        gpsStart = time[0]
        gps_times, cleaned, sample_mask, fs_new, breakdown = process_band(
            strain,
            channel_dict,
            gpsStart,
            f_band,
            fap=fap,
            nperseg=nperseg,
            bound=bound,
            plot=False,
            loss_stats=True,
        )
        save_hdf5_file(
            process_dirname,
            files.name,
            f_band,
            gps_times,
            cleaned,
            sample_mask,
            fs_new,
            breakdown,
        )
