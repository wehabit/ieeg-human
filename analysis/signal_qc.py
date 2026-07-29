"""Shared signal-quality helpers used by the active LC-proxy pipelines."""

import numpy as np
from scipy import ndimage, signal


def robust_z(v):
    med = np.median(v); mad = np.median(np.abs(v - med)) + 1e-12
    return (v - med) / (1.4826 * mad)


def dilate_boolean_mask(mask, half_width):
    """Boolean equivalent of convolution with an odd all-ones window.

    ``np.convolve`` is quadratic in the window width for this use and made a five-second dilation
    at 2048 Hz dominate the multi-contact cache runtime.  A one-dimensional maximum filter
    performs the same centred inclusive dilation in linear time.
    """
    mask = np.asarray(mask, bool)
    half_width = int(half_width)
    if half_width <= 0 or not mask.any():
        return mask.copy()
    return ndimage.maximum_filter1d(
        mask.astype(np.uint8), size=2 * half_width + 1,
        mode="constant", cval=0).astype(bool)


def ied_clean_mask(y, sf, z_hf=5.0, z_amp=8.0, pad_s=0.5):
    """True where the channel is FREE of interictal epileptiform discharges / sharp artifact.
    Detected on 20-80 Hz (spikes are broadband-sharp) so the 11-16 Hz sigma band itself is not
    used as the detector; plus a raw-amplitude guard. Flags are dilated by +/- pad_s."""
    sos = signal.butter(4, [20.0, 80.0], btype="band", fs=sf, output="sos")
    hf = np.abs(signal.hilbert(signal.sosfiltfilt(sos, y)))
    bad = (robust_z(hf) > z_hf) | (np.abs(robust_z(y)) > z_amp)
    k = int(pad_s * sf)
    bad = dilate_boolean_mask(bad, k)
    return ~bad


__all__ = ["robust_z", "dilate_boolean_mask", "ied_clean_mask"]
