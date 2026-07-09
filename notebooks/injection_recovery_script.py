from __future__ import annotations

import attrs
import numpy as np
from pycbc.detector import Detector
from tqdm.auto import tqdm


@attrs.define(kw_only=True)
class InjectionRecovery:
    fs: int
    f_band: tuple[float, float]
    f_sig: float
    amplitude: float
    sky_pos: tuple[float, float]
    gps_start: float
    duration: float
    detector_name: str
    det: Detector = attrs.field(init=False)

    def __attrs_post_init__(self):
        self.det = Detector(self.detector_name)

    @property
    def time_array(self) -> np.ndarray:
        return self.gps_start + np.arange(self.fs * self.duration) / self.fs

    @property
    def time_array_ref(self) -> np.ndarray:
        t_ref = self.time_array[len(self.time_array) // 2]
        return self.time_array - t_ref

    def antenna_pattern(self) -> tuple[np.ndarray, np.ndarray]:
        return self.det.antenna_pattern(
            self.sky_pos[0], self.sky_pos[1], 0, self.time_array
        )  # pyright: ignore[reportReturnType]

    def inject_signal(self) -> np.ndarray:
        fplus, fcross = self.antenna_pattern()
        antenna_amplitude = fplus + 1j * fcross
        return (
            self.amplitude
            * antenna_amplitude
            * np.exp(1j * (2 * np.pi * self.f_sig * self.time_array_ref))
        )

    def generate_noise(self) -> np.ndarray:
        return np.random.normal(
            0, np.sqrt(0.5), len(self.time_array)
        ) + 1j * np.random.normal(0, np.sqrt(0.5), len(self.time_array))

    def covariance_mat(self) -> np.ndarray:
        fplus, fcross = self.antenna_pattern()
        vp = np.dot(fplus, fplus)
        vc = np.dot(fcross, fcross)
        vpc = np.dot(fplus, fcross)
        return np.array([[vp, vpc], [vpc, vc]])

    def data_projection_fourier(self, data: np.ndarray) -> np.ndarray:
        dplus = np.fft.fft(data * self.antenna_pattern()[0])
        dcross = np.fft.fft(data * self.antenna_pattern()[1])
        return np.stack([dplus, dcross], axis=0)

    def scores(self, data: np.ndarray) -> np.ndarray:
        return np.einsum(
            "in,ij,jn->n",
            self.data_projection_fourier(data).conj(),
            np.linalg.inv(self.covariance_mat()),
            self.data_projection_fourier(data),
        )

    @property
    def optimal_snr_sq(self) -> float:
        return np.sum(np.abs(self.inject_signal()) ** 2)

    def recovered_scores(self, n_sims: int) -> np.ndarray:
      f_idx = int(round(self.f_sig * self.duration))
      scores = np.zeros(n_sims, dtype=np.float64)
      for i in tqdm(range(n_sims), desc="Running simulations"):
          noise = self.generate_noise()
          data = self.inject_signal() + noise
          scores[i] = self.scores(data)[f_idx].real
      return scores
