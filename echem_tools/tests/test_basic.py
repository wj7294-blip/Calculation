import sys
from pathlib import Path
import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from echem.config import load_config
from echem.eis.synth import synth_spectrum, circuit_truth
from echem.eis.drt import compute_drt_builtin
from echem.eis.preprocess import preprocess_eis, linkk_check
from echem.io.loaders import load_eis, load_chrono
from echem.floating.synth import _cell_curve, write_mpt, write_dta, write_neware_csv
from echem.floating.analysis import analyze_cell
from echem.pipeline_eis import _write_synthetic_eis


@pytest.fixture(scope="module")
def cfg():
    return load_config(ROOT / "config.yaml")


def test_config_floats(cfg):
    assert isinstance(cfg["eis"]["freq_max_Hz"], float)


def test_drt_recovers_two_rc(cfg):
    spec = {"R0": 5.0, "RC": [[10.0, 1e-4], [20.0, 1e-2]], "L": 1e-7}
    d = synth_spectrum(spec, 1e5, 1e-2, 10, 1.0, np.random.default_rng(3))
    d, _ = preprocess_eis(d, cfg["eis"])
    r = compute_drt_builtin(d, cfg["eis"]["drt"])
    assert abs(r.R_inf - 5.0) < 0.2
    assert len(r.peaks) == 2
    for (R, tau), (_, p) in zip(spec["RC"], r.peaks.iterrows()):
        assert abs(np.log10(p.tau_s / tau)) < 0.1
        assert abs(p.R_ohm / R - 1) < 0.08


def test_manual_lambda_merges_close_peaks(cfg):
    spec = {"R0": 4.0, "RC": [[10.0, 2e-3], [10.0, 1.6e-2]]}
    d = synth_spectrum(spec, 1e5, 1e-2, 10, 1.0, np.random.default_rng(1))
    r_small = compute_drt_builtin(d, cfg["eis"]["drt"], lam=1e-5)
    r_big = compute_drt_builtin(d, cfg["eis"]["drt"], lam=1.0)
    assert len(r_small.peaks) == 2 and len(r_big.peaks) == 1


def test_linkk_clean_vs_noisy(cfg):
    spec = {"R0": 5.0, "RC": [[10.0, 1e-4], [20.0, 1e-2]]}
    d = synth_spectrum(spec, 1e5, 1e-2, 10, 0.2, np.random.default_rng(2))
    kk = linkk_check(d, cfg["eis"]["kk"])
    assert kk.rms_res_pct < 1.0 and not kk.warn
    bad = d
    bad.Z = bad.Z * (1 + 0.05 * np.linspace(0, 1, len(bad.Z)))   # drift -> KK violation
    kk2 = linkk_check(bad, cfg["eis"]["kk"])
    assert kk2.rms_res_pct > kk.rms_res_pct


def test_eis_loaders_roundtrip(tmp_path, cfg):
    spec = {"R0": 5.0, "RC": [[10.0, 1e-4]]}
    d = synth_spectrum(spec, 1e5, 1e-1, 8, 0.0, np.random.default_rng(0))
    for i in range(3):
        fn = _write_synthetic_eis(d, tmp_path, f"s{i}", i)
        d2 = load_eis(tmp_path / fn, cfg, "generic" if fn.endswith(".csv") else None)
        assert len(d2.freq) == len(d.freq)
        assert np.allclose(d2.Z.imag, d.Z.imag, rtol=1e-4, atol=1e-5)


def test_chrono_loaders_and_metrics(tmp_path, cfg):
    rng = np.random.default_rng(0)
    mass_g = 0.008
    df = _cell_curve(rng, 0.8 * mass_g, mass_g, cfg["floating"], cfg["synthetic"]["floating"])
    write_mpt(df, tmp_path / "a.mpt"); write_dta(df, tmp_path / "a.DTA"); write_neware_csv(df, tmp_path / "a.csv")
    vals = []
    for fn, prof in [("a.mpt", None), ("a.DTA", None), ("a.csv", "neware")]:
        cd = load_chrono(tmp_path / fn, cfg, prof)
        r = analyze_cell(cd, fn, 8.0, 1.13, cfg["floating"])
        vals.append(r.metrics["leak_mA_g"])
        assert 23 < r.metrics["cv_duration_h"] < 25
        assert r.metrics["Q_leak_ex_transient_mAh_g"] < r.metrics["Q_leak_mAh_g"]
        assert r.metrics["decay_slope"] < 0
    assert np.allclose(vals, vals[0], rtol=1e-3)
    assert abs(vals[0] / 0.8 - 1) < 0.3


def test_two_way_anova_detects_interaction():
    from echem.floating.stats import group_stats
    rng = np.random.default_rng(0)
    rows = []
    for s in ["bare", "SAM"]:
        for e in ["A", "B"]:
            mu = 1.0 - 0.3 * (s == "SAM") + 0.6 * (s == "SAM" and e == "B")   # crossover interaction
            for i in range(6):
                rows.append({"surface": s, "electrolyte": e, "leak_mA_g": mu + 0.05 * rng.standard_normal()})
    cells = pd.DataFrame(rows)
    st = group_stats(cells, {"factors": ["surface", "electrolyte"], "metrics": ["leak_mA_g"], "nonparametric": "always"})
    a = st["metrics"]["leak_mA_g"]["anova"]
    assert a.loc["surface:electrolyte", "PR(>F)"] < 0.01
    assert st["metrics"]["leak_mA_g"]["nonparametric"]["srh"].loc["surface:electrolyte", "p"] < 0.05
