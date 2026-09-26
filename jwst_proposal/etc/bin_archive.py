"""Estimate R~100 S/N from an ETC archive, assuming diagonal spectral covariance.

Run: python jwst_proposal/etc/bin_archive.py ARCHIVE OUTPUT_DIRECTORY
This reads files directly from the tar; it does not extract arbitrary paths.
It does not run Pandeia or predict calibration/host-decomposition accuracy.
"""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np
from astropy.io import fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--centers", nargs="+", type=float,
                        help="Observed bin centers in microns; defaults to reference wavelength")
    args = parser.parse_args()
    with tarfile.open(args.archive) as archive:
        def read(suffix):
            matches = [m for m in archive.getmembers()
                       if m.isfile() and m.name.endswith("/" + suffix)]
            if len(matches) != 1:
                raise ValueError(f"Expected exactly one {suffix}")
            return archive.extractfile(matches[0]).read()

        raw_input = read("input.json")
        config = json.loads(raw_input)

        def spectrum(name):
            with fits.open(io.BytesIO(read(f"lineplot/lineplot_{name}.fits"))) as hdus:
                return (np.array(hdus[1].data.field(0)),
                        np.array(hdus[1].data.field(1)))

        wave, flux = spectrum("extracted_flux")
        noise_wave, noise = spectrum("extracted_noise")
        sn_wave, sn = spectrum("sn")
        if not (np.array_equal(wave, noise_wave) and np.array_equal(wave, sn_wave)):
            raise ValueError("Wavelength grids differ")
        if not (np.all(np.diff(wave) > 0) and np.all(noise > 0)
                and np.all(np.isfinite(flux)) and np.all(np.isfinite(noise))):
            raise ValueError("Invalid spectral values")
        if not np.allclose(flux / noise, sn, rtol=1e-8, atol=1e-10):
            raise ValueError("Flux/noise does not reproduce ETC S/N")

    # Select whole output pixels whose centers lie in lambda +/- lambda/(2R).
    # No interpolation into extra samples and no claim of a supplied covariance.
    pixel_edges = np.r_[wave[0] - (wave[1] - wave[0]) / 2,
                        (wave[1:] + wave[:-1]) / 2,
                        wave[-1] + (wave[-1] - wave[-2]) / 2]
    def bin_at(center):
        lo, hi = center * (1 - 1 / 200), center * (1 + 1 / 200)
        if lo < pixel_edges[0] or hi > pixel_edges[-1]:
            return None
        idx = np.flatnonzero((wave >= lo) & (wave < hi))
        signal = float(flux[idx].sum())
        sigma = float(np.sqrt(np.square(noise[idx]).sum()))
        width = float(pixel_edges[idx[-1] + 1] - pixel_edges[idx[0]])
        return dict(center_um=center, requested_R=100,
                         actual_R=center / width, n_pixels=len(idx),
                         summed_signal_e_s=signal, diagonal_noise_e_s=sigma,
                         snr_diagonal=signal / sigma)

    anchor = config["strategy"]["reference_wavelength"]
    centers = args.centers if args.centers is not None else [anchor]
    rows = [row for center in centers
            if (row := bin_at(center)) is not None]
    # Contiguous nominal R=100 bins anchored on the reference wavelength.
    # Shared boundaries assign each whole pixel to exactly one bin.
    ratio = 201 / 199
    low_index = int(np.floor(np.log(wave[0] / anchor) / np.log(ratio))) - 1
    high_index = int(np.ceil(np.log(wave[-1] / anchor) / np.log(ratio))) + 1
    full_rows = [row for k in range(low_index, high_index + 1)
                 if (row := bin_at(anchor * ratio**k)) is not None]
    if not rows or not full_rows:
        raise ValueError("No complete R=100 bins at requested centers or in sub-band")
    summary = dict(
        archive=args.archive.name,
        archive_sha256=hashlib.sha256(args.archive.read_bytes()).hexdigest(),
        configuration=config["configuration"],
        scene=config["scene"], strategy=config["strategy"],
        wavelength_range_um=[float(wave[0]), float(wave[-1])],
        output_pixel_count=len(wave),
        native_snr_at_reference=float(np.interp(
            config["strategy"]["reference_wavelength"], wave, sn)),
        method="Whole-pixel sum(flux)/sqrt(sum(noise**2)); diagonal covariance assumption",
        limitations=["No spectral covariance supplied or propagated",
                     "No added calibration, resampling, or host-decomposition errors",
                     "One point source only; host photon noise is absent",
                     "No engine version, sky coordinates, or scalar timing report in input.json",
                     "One channel/sub-band only; does not establish full program time"],
        bins=rows,
        full_subband_bins=full_rows,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "input.json").write_bytes(raw_input)
    (args.output / "analysis.json").write_text(json.dumps(summary, indent=2) + "\n")
    with (args.output / "r100_estimates.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (args.output / "r100_full_subband.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(full_rows[0]))
        writer.writeheader()
        writer.writerows(full_rows)
    print(f"Full sub-band: {len(full_rows)} complete bins; diagonal S/N "
          f"{min(r['snr_diagonal'] for r in full_rows):.2f} to "
          f"{max(r['snr_diagonal'] for r in full_rows):.2f}")
    for row in rows:
        print(f"{row['center_um']:.2f} um: approximate R=100 S/N "
              f"{row['snr_diagonal']:.2f}, {row['n_pixels']} pixels")


if __name__ == "__main__":
    main()
