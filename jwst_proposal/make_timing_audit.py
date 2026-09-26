"""Summarize saved histories; do not infer event onsets or reverberation lags."""
import json
from pathlib import Path

import numpy as np
from astropy.time import Time

HERE = Path(__file__).resolve().parent
INPUTS = HERE / 'inputs'


def build_audit():
    snapshot = json.loads((INPUTS / 'timing_history_snapshot.json').read_text())
    targets = []
    for target in snapshot['targets']:
        bands = {}
        for name, series in target['bands'].items():
            data = np.asarray(series['data'], dtype=float)
            data = data[np.isfinite(data[:, :2]).all(axis=1)]
            years = Time(data[:, 0], format='mjd').decimalyear
            bands[name] = dict(
                units=series['units'], points=len(data),
                start_utc=Time(data[:, 0].min(), format='mjd').isot[:10],
                end_utc=Time(data[:, 0].max(), format='mjd').isot[:10],
                annual_medians=[dict(year=int(year), points=int(np.sum(years.astype(int) == year)),
                                     median=float(np.median(data[years.astype(int) == year, 1])))
                                for year in sorted(set(years.astype(int)))],
            )
        targets.append(dict(name=target['name'], z=target['z'], bands=bands,
                            spectral_dates=target['spectral_dates'],
                            existing_review=target['existing_review'],
                            spectral_review=target['spectral_review'],
                            extraction=target.get('extraction', {}),
                            expansion_rationale=target.get('expansion_rationale')))
    return dict(
        description='Timing evidence from saved light curves; calendar-year medians without inter-survey rescaling or host subtraction. Display-series counts are not raw exposure counts.',
        cycle6_nominal_dates=['2027-07-01', '2028-06-30'],
        event_onsets_fitted=False, target_dust_lags_measured=False,
        normalized_response_phase_coverage_validated=False,
        interpretation='Years covered by photometry or spectroscopy do not themselves determine time since an illumination change. Infer event intervals and luminosity-dependent lag distributions before asserting response-phase coverage.',
        selection_decisions={
            'P9694': 'Retain as a variable science candidate: optical rise between the 2021 and 2022 annual medians, brighter through 2025, W1 rising through 2024. Not a stable control.',
            'P9584': 'Remove from the working request. Small optical changes and no demonstrated matched-control role; original data remain in this decision audit.',
        },
        observed_trend_intervals={
            'P1823': 'Sustained optical rise in 2019–2023; W1 minimum around 2018–2019 and recovery through 2024; spectral high state in September 2026.',
            'P11530': 'Optical decline in 2018–2021, partial recovery in 2022–2023, renewed fading through 2025. W1 also shows decline and recovery.',
            'P9694': 'Optical rise between 2021 and 2022 annual medians; brighter through 2025, while W1 rises through 2024.',
            'P11113': 'Optical r brightening between 2021 and 2022, decline through 2025; W1 rises through 2024.',
        },
        pending=['Bracket optical change times and amplitudes with calibrated, host-corrected data; spectral observation dates alone are not transition dates.',
                 'Estimate rest-frame dust lag distributions and luminosity uncertainty; assess event-age/lag coverage at allowed JWST dates.',
                 'Check historical optical coverage of the three literature anchors; two do not have full optical display series in this snapshot.',
                 'Obtain contemporary optical coverage to bridge saved light curves into the JWST epoch.',
                 'Validate the inter-instrument separation and any repeat-visit benefit with source-specific predictions and APT visibility.'],
        targets=targets,
    )


if __name__ == '__main__':
    audit = build_audit()
    (INPUTS / 'timing_audit.json').write_text(json.dumps(audit, indent=2)+'\n')
    print(f"Timing audit: {len(audit['targets'])} histories; response-phase coverage remains unvalidated.")
