"""Run the project-local MIRI/MRS ETC and optionally compare a web export.

Use .venv/bin/python run_local.py INPUT_JSON OUTPUT_DIR [--compare-web TAR].
Reference paths are resolved relative to this script; no shell setup is needed.
"""

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import time


def configure():
    base = Path(__file__).resolve().parent
    paths = {
        "pandeia_refdata": base / "reference_data/pandeia_data-2026.7-jwst",
        "PSF_DIR": base / "reference_data/mrs_psfs_2026.7",
        "PYSYN_CDBS": base / "reference_data/synphot",
    }
    for key, path in paths.items():
        if not path.is_dir():
            raise FileNotFoundError(f"Missing {key}: {path}")
        os.environ[key] = str(path)
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--compare-web", type=Path)
    args = parser.parse_args()
    paths = configure()
    import numpy as np
    from astropy.io import fits
    from synphot.config import conf as synphot_conf
    synphot_conf.vega_file = str(paths["PYSYN_CDBS"] / "calspec/alpha_lyr_stis_010.fits")
    from pandeia.engine import __version__, pandeia_version
    from pandeia.engine.perform_calculation import perform_calculation

    def plain(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(type(value).__name__)

    raw = args.input.read_bytes()
    config = json.loads(raw)
    instrument = config["configuration"]["instrument"]
    if instrument["instrument"] != "miri" or instrument["mode"] != "mrs":
        raise ValueError("This installation contains MIRI/MRS PSFs only")
    if args.output.exists():
        raise FileExistsError(f"Choose a new output directory: {args.output}")
    pandeia_version()
    started = time.monotonic()
    report = perform_calculation(config, dict_report=False)
    result = report.as_dict()
    elapsed = time.monotonic() - started
    summary = {
        "engine_version": __version__,
        "reference_data_version": (paths["pandeia_refdata"] / "VERSION_DATA").read_text().splitlines()[0],
        "psf_version": (paths["PSF_DIR"] / "VERSION_PSF").read_text().splitlines()[0],
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "runtime_seconds": elapsed,
        "scalar": result["scalar"], "warnings": result["warnings"],
    }
    if args.compare_web:
        comparisons = {}
        with tarfile.open(args.compare_web) as archive:
            def member(suffix):
                found = [m for m in archive.getmembers()
                         if m.isfile() and m.name.endswith("/" + suffix)]
                if len(found) != 1:
                    raise ValueError(f"Expected one {suffix} in web archive")
                return archive.extractfile(found[0]).read()

            if json.loads(member("input.json")) != config:
                raise ValueError("Input differs from the web calculation")
            for key in ("extracted_flux", "extracted_noise", "sn"):
                with fits.open(io.BytesIO(member(f"lineplot/lineplot_{key}.fits"))) as hdus:
                    wave = np.array(hdus[1].data.field(0))
                    expected = np.array(hdus[1].data.field(1))
                actual_wave, actual = result["1d"][key]
                np.testing.assert_allclose(actual_wave, wave, rtol=1e-10, atol=1e-12)
                relative = np.abs(actual - expected) / np.maximum(np.abs(expected), 1e-30)
                comparisons[key] = {
                    "samples": len(wave),
                    "max_relative_difference": float(np.max(relative)),
                    "matches_rtol_1e_5": bool(np.allclose(actual, expected, rtol=1e-5, atol=1e-10)),
                }
        summary["web_validation"] = {
            "archive": args.compare_web.name,
            "archive_sha256": hashlib.sha256(args.compare_web.read_bytes()).hexdigest(),
            "input_identical": True, "spectra": comparisons,
            "passed": all(item["matches_rtol_1e_5"] for item in comparisons.values()),
        }

    args.output.mkdir(parents=True)
    (args.output / "input.json").write_bytes(raw)
    (args.output / "report.json").write_text(json.dumps(summary, default=plain, indent=2) + "\n")
    line_dir = args.output / "lineplot"
    line_dir.mkdir()
    for key, hdus in report.as_fits()["1d"].items():
        if not isinstance(hdus[0], fits.PrimaryHDU):
            hdus = fits.HDUList([fits.PrimaryHDU(), *list(hdus)])
        hdus.writeto(line_dir / f"lineplot_{key}.fits")
    with tarfile.open(args.output / "local_export.tar", "w") as archive:
        archive.add(args.output / "input.json", arcname="local/input.json")
        archive.add(line_dir, arcname="local/lineplot")
    print(json.dumps(summary, default=plain, indent=2))
    if summary.get("web_validation", {}).get("passed") is False:
        raise RuntimeError("Local spectra differ from the web export; inspect report.json")


if __name__ == "__main__":
    main()
