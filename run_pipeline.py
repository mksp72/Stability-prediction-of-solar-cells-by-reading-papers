"""Run the existing pipeline.py on real CSV data (no synthetic data)."""
import argparse
from pathlib import Path
import pandas as pd
from pipeline import validate, fit_post_peak


def main():
    parser = argparse.ArgumentParser(description='Validate and fit perovskite stability curves.')
    parser.add_argument('--devices', default='devices.csv', help='Device metadata CSV')
    parser.add_argument('--curves', default='curves.csv', help='Time-retention CSV')
    parser.add_argument('--output', default='results', help='Output folder')
    args = parser.parse_args()
    device_path, curve_path = Path(args.devices), Path(args.curves)
    missing = [str(p) for p in (device_path, curve_path) if not p.is_file()]
    if missing:
        print('Missing input files: ' + ', '.join(missing))
        print('Required devices.csv columns: device_id, doi, pce0, voc0, jsc0, ff0')
        print('Required curves.csv columns: device_id, protocol, time_h, retention_pct')
        print('Place both CSV files next to this runner or pass --devices and --curves with full paths.')
        raise SystemExit(1)
    devices = pd.read_csv(device_path)
    curves = pd.read_csv(curve_path)
    validate(devices, curves)
    print(f'Validated {len(devices)} devices and {len(curves)} curve data points.')
    rows = []
    for (device_id, protocol), part in curves.groupby(['device_id', 'protocol'], sort=True):
        if protocol != 'ISOS-L':
            continue
        for n in (2, 3):
            item = {'device_id': device_id, 'protocol': protocol, 'n_components': n}
            try:
                fit = fit_post_peak(part.time_h.to_numpy(), part.retention_pct.to_numpy(), n_components=n)
                item.update({'status': 'ok', 'peak_time_h': fit.peak_time_h,
                    'peak_retention_pct': fit.peak_retention_pct,
                    'observed_t80_h': fit.observed_t80_h,
                    'fitted_t80_h': fit.fitted_t80_h,
                    'fitted_t80_extrapolated': fit.fitted_t80_extrapolated,
                    'rmse_post_peak_pct': fit.rmse_post_peak_pct})
                for j, (amp, rate) in enumerate(zip(fit.amplitudes, fit.rates_per_h), 1):
                    item[f'amplitude_{j}'] = amp
                    item[f'rate_{j}_per_h'] = rate
            except (ValueError, ZeroDivisionError, FloatingPointError) as exc:
                item.update({'status': f'not_fitted: {exc}'})
            rows.append(item)
    if not rows:
        print('No ISOS-L curves found: no fitting output created.')
        return
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    destination = output / 'isos_l_fit_results.csv'
    pd.DataFrame(rows).to_csv(destination, index=False)
    print(f'Wrote {len(rows)} fitting rows to {destination.resolve()}')
    print('Note: these fits are exploratory; triple-exponential parameters may not be identifiable.')


if __name__ == '__main__':
    main()
