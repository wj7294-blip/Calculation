"""File loaders: BioLogic .mpt, Gamry .DTA, generic/Neware/LANHE CSV.

Every loader returns either
  * EISData   : freq [Hz], Z [ohm, complex, Z = Z' + jZ'' with Z'' < 0 for capacitive], extra per-point columns
  * ChronoData: DataFrame with standard columns  time_s, voltage_V, current_mA  (+ step, cycle if available)
Current sign convention: positive = charging (oxidation of the working electrode).
"""
from __future__ import annotations
import csv
import io
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


class LoaderError(RuntimeError):
    pass


@dataclass
class EISData:
    freq: np.ndarray
    Z: np.ndarray
    source: str = ""
    extra: pd.DataFrame | None = None     # time, cycle, Ewe ... per point
    meta: dict = field(default_factory=dict)

    def sort(self) -> "EISData":
        idx = np.argsort(self.freq)[::-1]
        ex = self.extra.iloc[idx].reset_index(drop=True) if self.extra is not None else None
        return EISData(self.freq[idx], self.Z[idx], self.source, ex, dict(self.meta))

    def to_frame(self) -> pd.DataFrame:
        df = pd.DataFrame({"freq_Hz": self.freq, "Zre_ohm": self.Z.real, "Zim_ohm": self.Z.imag})
        if self.extra is not None:
            for c in self.extra.columns:
                if c not in df:
                    df[c] = self.extra[c].values
        return df


@dataclass
class ChronoData:
    df: pd.DataFrame          # time_s, voltage_V, current_mA [, step, cycle]
    source: str = ""
    meta: dict = field(default_factory=dict)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def _read_text(path: Path) -> str:
    raw = path.read_bytes()
    for enc in ("utf-8-sig", "latin-1", "cp1252"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("latin-1", errors="replace")


def _to_float_frame(df: pd.DataFrame) -> pd.DataFrame:
    """Convert numeric-looking columns (handles decimal comma)."""
    out = {}
    for c in df.columns:
        s = df[c]
        if not pd.api.types.is_numeric_dtype(s):
            s2 = pd.to_numeric(s.astype(str).str.replace(",", ".", regex=False).str.strip(), errors="coerce")
            if s2.notna().mean() > 0.8:
                s = s2
        out[c] = s
    return pd.DataFrame(out)


def _find_col(columns, candidates) -> Optional[str]:
    norm = {re.sub(r"\s+", "", c).lower(): c for c in columns}
    for cand in candidates:
        key = re.sub(r"\s+", "", cand).lower()
        if key in norm:
            return norm[key]
    # fuzzy: startswith
    for cand in candidates:
        key = re.sub(r"\s+", "", cand).lower()
        for k, c in norm.items():
            if k.startswith(key):
                return c
    return None


def detect_kind(df: pd.DataFrame) -> str:
    """'eis' if impedance-like columns present, else 'chrono'."""
    cols = [c.lower() for c in df.columns]
    if any(("freq" in c) for c in cols) and any(("re(z)" in c or "zreal" in c or "z'" in c or "zre" in c) for c in cols):
        return "eis"
    return "chrono"


# ----------------------------------------------------------------------------
# BioLogic .mpt
# ----------------------------------------------------------------------------
MPT_EIS_MAP = {
    "freq_Hz": ["freq/Hz"],
    "Zre_ohm": ["Re(Z)/Ohm"],
    "negZim_ohm": ["-Im(Z)/Ohm"],
    "time_s": ["time/s"],
    "voltage_V": ["<Ewe>/V", "Ewe/V"],
    "cycle": ["cycle number"],
}
MPT_CHRONO_MAP = {
    "time_s": ["time/s"],
    "voltage_V": ["Ewe/V", "<Ewe>/V", "Ecell/V"],
    "current_mA": ["<I>/mA", "I/mA", "control/mA"],
    "step": ["Ns", "mode"],
    "cycle": ["cycle number", "half cycle"],
}


def read_mpt(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    text = _read_text(path)
    lines = text.splitlines()
    n_header = None
    for ln in lines[:5]:
        m = re.search(r"Nb header lines\s*:\s*(\d+)", ln)
        if m:
            n_header = int(m.group(1))
            break
    if n_header is None:
        # fall back: first line containing 'freq/Hz' or 'time/s'
        for i, ln in enumerate(lines):
            if "\t" in ln and ("time/s" in ln or "freq/Hz" in ln):
                n_header = i + 1
                break
    if n_header is None:
        raise LoaderError(f"{path}: cannot find header line count")
    header = lines[n_header - 1].rstrip("\n").split("\t")
    header = [h.strip() for h in header]
    body = "\n".join(lines[n_header:])
    df = pd.read_csv(io.StringIO(body), sep="\t", header=None, names=header, engine="python",
                     usecols=range(len(header)), index_col=False)
    df = _to_float_frame(df)
    df.attrs["header_lines"] = lines[:n_header]
    return df


# ----------------------------------------------------------------------------
# Gamry .DTA
# ----------------------------------------------------------------------------
def read_dta(path: str | Path) -> dict[str, pd.DataFrame]:
    """Return dict of tables {name: DataFrame}. Tables: ZCURVE, CURVE, CURVE1, OCVCURVE ..."""
    path = Path(path)
    lines = _read_text(path).splitlines()
    tables: dict[str, pd.DataFrame] = {}
    i = 0
    while i < len(lines):
        parts = lines[i].split("\t")
        if len(parts) >= 2 and parts[1].strip().upper() == "TABLE":
            name = parts[0].strip()
            header = [h.strip() for h in lines[i + 1].split("\t")]
            # units line follows header
            j = i + 3
            rows = []
            while j < len(lines):
                row = lines[j].split("\t")
                if len(row) < 2 or (len(row) >= 2 and row[1].strip().upper() == "TABLE"):
                    break
                # data rows start with a tab (empty first field) in Gamry files
                rows.append(row)
                j += 1
            if rows:
                width = max(len(header), max(len(r) for r in rows))
                header_full = header + [f"col{k}" for k in range(len(header), width)]
                df = pd.DataFrame([r + [""] * (width - len(r)) for r in rows], columns=header_full)
                # drop leading empty column (Gamry prefixes each data row with a tab)
                if df.columns[0] == "" or (df.iloc[:, 0].astype(str).str.strip() == "").all():
                    df = df.iloc[:, 1:]
                    if header[0] == "":
                        pass
                df = _to_float_frame(df)
                # Gamry numbers columns Pt, Time, ... ; after dropping blank col headers may be shifted
                tables[name] = df.reset_index(drop=True)
            i = j
        else:
            i += 1
    if not tables:
        raise LoaderError(f"{path}: no TABLE blocks found")
    return tables


def _dta_fix_columns(df: pd.DataFrame, expected: list[str]) -> pd.DataFrame:
    """Gamry header row starts with the column names, and data rows start with an empty
    field, which shifts names by one. Align by expected first names."""
    cols = list(df.columns)
    if "Pt" in cols and cols.index("Pt") != 0:
        # header had a leading blank: columns are already aligned after dropping blank col
        pass
    if cols and cols[0] == "" and len(cols) > 1:
        df = df.iloc[:, 1:]
    return df


# ----------------------------------------------------------------------------
# generic CSV with profile mapping
# ----------------------------------------------------------------------------
def read_csv_profile(path: str | Path, profile: dict) -> pd.DataFrame:
    path = Path(path)
    text = _read_text(path)
    delim = profile.get("delimiter")
    # skip leading comment/metadata lines until a line with >=2 delimiters
    lines = text.splitlines()
    if delim is None:
        sample = "\n".join(lines[:50])
        try:
            delim = csv.Sniffer().sniff(sample, delimiters=",;\t").delimiter
        except csv.Error:
            delim = ","
    start = 0
    for i, ln in enumerate(lines):
        if ln.count(delim) >= 2 and not ln.startswith("#"):
            start = i
            break
    df = pd.read_csv(io.StringIO("\n".join(lines[start:])), sep=delim, engine="python", comment="#")
    df.columns = [str(c).strip() for c in df.columns]
    return _to_float_frame(df)


def map_columns(df: pd.DataFrame, profile: dict) -> pd.DataFrame:
    """Map raw columns -> standard names using profile {std: {candidates, scale}}."""
    out = pd.DataFrame(index=df.index)
    found = {}
    for std, spec in profile.items():
        if not isinstance(spec, dict) or "candidates" not in spec:
            continue
        c = _find_col(df.columns, spec["candidates"])
        if c is not None:
            scale = float(spec.get("scale", 1.0))
            s = df[c]
            out[std] = s * scale if pd.api.types.is_numeric_dtype(s) else s
            found[std] = c
    if "current_mA" not in out and "current_A" in out:
        out["current_mA"] = out["current_A"]
    if "Zim_ohm" not in out and "negZim_ohm" in out:
        out["Zim_ohm"] = -out["negZim_ohm"]
    out.attrs["column_map"] = found
    return out


# ----------------------------------------------------------------------------
# public API
# ----------------------------------------------------------------------------
def load_eis(path: str | Path, cfg: dict | None = None, profile_name: str | None = None) -> EISData:
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".mpt":
        raw = read_mpt(path)
        std = map_columns(raw, {k: {"candidates": v} for k, v in MPT_EIS_MAP.items()})
        meta = {"format": "biologic_mpt"}
    elif ext == ".dta":
        tables = read_dta(path)
        key = next((k for k in tables if k.upper().startswith("ZCURVE")), None)
        if key is None:
            raise LoaderError(f"{path}: no ZCURVE table (not an EIS file?)")
        raw = tables[key]
        std = map_columns(raw, {
            "freq_Hz": {"candidates": ["Freq"]}, "Zre_ohm": {"candidates": ["Zreal"]},
            "Zim_ohm": {"candidates": ["Zimag"]}, "time_s": {"candidates": ["Time", "T"]},
            "voltage_V": {"candidates": ["Vdc"]},
        })
        meta = {"format": "gamry_dta"}
    else:
        raw, std, used = _csv_auto(path, cfg, profile_name, ("freq_Hz", "Zre_ohm", "Zim_ohm"))
        meta = {"format": f"csv:{used}"}
    for need in ("freq_Hz", "Zre_ohm", "Zim_ohm"):
        if need not in std:
            raise LoaderError(f"{path}: column for {need} not found (have {list(raw.columns)[:12]}...)")
    std = std.dropna(subset=["freq_Hz", "Zre_ohm", "Zim_ohm"])
    std = std[std["freq_Hz"] > 0]
    f = std["freq_Hz"].to_numpy(float)
    Z = std["Zre_ohm"].to_numpy(float) + 1j * std["Zim_ohm"].to_numpy(float)
    extra = std.drop(columns=[c for c in ["freq_Hz", "Zre_ohm", "Zim_ohm", "negZim_ohm"] if c in std]).reset_index(drop=True)
    meta["column_map"] = std.attrs.get("column_map", {})
    return EISData(f, Z, str(path), extra if len(extra.columns) else None, meta).sort()


def load_chrono(path: str | Path, cfg: dict | None = None, profile_name: str | None = None) -> ChronoData:
    path = Path(path)
    ext = path.suffix.lower()
    meta: dict = {}
    if ext == ".mpt":
        raw = read_mpt(path)
        std = map_columns(raw, {k: {"candidates": v} for k, v in MPT_CHRONO_MAP.items()})
        meta["format"] = "biologic_mpt"
    elif ext == ".dta":
        tables = read_dta(path)
        keys = [k for k in tables if k.upper().startswith("CURVE") or k.upper().startswith("OCVCURVE")]
        if not keys:
            raise LoaderError(f"{path}: no CURVE table")
        frames = []
        for k in keys:
            t = tables[k].copy()
            t["step"] = k
            frames.append(t)
        raw = pd.concat(frames, ignore_index=True)
        std = map_columns(raw, {
            "time_s": {"candidates": ["T", "Time"]}, "voltage_V": {"candidates": ["Vf", "V"]},
            "current_mA": {"candidates": ["Im", "I"], "scale": 1000.0},   # Gamry current in A
            "step": {"candidates": ["step"]},
        })
        meta["format"] = "gamry_dta"
    else:
        raw, std, used = _csv_auto(path, cfg, profile_name, ("time_s", "voltage_V", "current_mA"))
        meta["format"] = f"csv:{used}"
        meta["cv_step_labels"] = _profile(cfg, used).get("cv_step_labels", [])
    for need in ("time_s", "voltage_V", "current_mA"):
        if need not in std:
            raise LoaderError(f"{path}: column for {need} not found (have {list(raw.columns)[:12]}...)")
    keep = [c for c in ["time_s", "voltage_V", "current_mA", "step", "cycle"] if c in std]
    df = std[keep].dropna(subset=["time_s", "voltage_V", "current_mA"]).reset_index(drop=True)
    # time may be relative per step in some exports -> enforce monotonic by cumulative fix
    t = df["time_s"].to_numpy(float)
    if len(t) > 1 and np.any(np.diff(t) < 0):
        d = np.diff(t)
        resets = np.where(d < 0)[0]
        offs = np.zeros_like(t)
        for r in resets:
            offs[r + 1:] += t[r] - t[r + 1] + (d[d > 0].min() if np.any(d > 0) else 1.0)
        df["time_s"] = t + offs
        meta["time_reset_fixed"] = int(len(resets))
    meta["column_map"] = std.attrs.get("column_map", {})
    return ChronoData(df, str(path), meta)


def _csv_auto(path: Path, cfg: dict | None, profile_name: str | None, required: tuple):
    """Read a CSV with the given profile, or - if none given - try the default profile first and
    then every other profile in config.yaml until the required standard columns are found."""
    cfg = cfg or {}
    if profile_name:
        order = [profile_name]
    else:
        default = cfg.get("default_csv_profile", "generic")
        order = [default] + [k for k in cfg.get("csv_profiles", {}) if k != default]
    last = None
    for name in order:
        prof = _profile(cfg, name)
        raw = read_csv_profile(path, prof)
        std = map_columns(raw, prof)
        if all(r in std for r in required):
            return raw, std, name
        last = (raw, std, name)
    return last if last else (pd.DataFrame(), pd.DataFrame(), order[0])


def _profile(cfg: dict | None, name: str | None) -> dict:
    cfg = cfg or {}
    profiles = cfg.get("csv_profiles", {})
    name = name or cfg.get("default_csv_profile", "generic")
    if name not in profiles:
        raise LoaderError(f"csv profile '{name}' not in config.yaml csv_profiles ({list(profiles)})")
    return profiles[name]


def load_file(path: str | Path, kind: str, cfg: dict | None = None, profile_name: str | None = None):
    if kind == "eis":
        return load_eis(path, cfg, profile_name)
    if kind in ("float", "chrono", "floating"):
        return load_chrono(path, cfg, profile_name)
    raise ValueError(f"unknown kind {kind}")
