from __future__ import annotations
from pathlib import Path
import pandas as pd

REQUIRED = ["cell_id", "file"]
STANDARD = ["cell_id", "file", "type", "surface", "electrolyte", "active_mass_mg", "area_cm2", "note"]


class MetadataError(ValueError):
    pass


def load_metadata(path: str | Path, data_dir: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise MetadataError(f"metadata file not found: {path}")
    df = pd.read_csv(path, comment="#", skip_blank_lines=True, dtype=str)
    df.columns = [c.strip() for c in df.columns]
    # Korean header alias
    df = df.rename(columns={"비고": "note", "메모": "note"})
    for c in REQUIRED:
        if c not in df.columns:
            raise MetadataError(f"metadata missing required column '{c}' (columns: {list(df.columns)})")
    df = df.dropna(subset=["cell_id", "file"])
    df = df[~df["cell_id"].astype(str).str.startswith("#")]
    for c in STANDARD:
        if c not in df.columns:
            df[c] = None
    for c in ["active_mass_mg", "area_cm2"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["file"] = df["file"].astype(str).str.strip()
    data_dir = Path(data_dir)
    df["path"] = [str((data_dir / f).resolve()) if not Path(f).is_absolute() else f for f in df["file"]]
    # infer type from folder if missing
    def _infer(row):
        t = row["type"]
        if isinstance(t, str) and t.strip():
            return t.strip().lower()
        p = row["file"].replace("\\", "/").lower()
        if "/eis" in p or p.startswith("eis"):
            return "eis"
        if "/float" in p or p.startswith("float"):
            return "float"
        return "unknown"
    df["type"] = df.apply(_infer, axis=1)
    missing = [p for p in df["path"] if not Path(p).exists()]
    df["exists"] = [Path(p).exists() for p in df["path"]]
    if missing:
        import warnings
        warnings.warn(f"{len(missing)} file(s) listed in metadata do not exist, e.g. {missing[:3]}")
    if df["cell_id"].duplicated().any():
        dup = df.loc[df["cell_id"].duplicated(), "cell_id"].tolist()
        raise MetadataError(f"duplicate cell_id in metadata: {dup[:5]}")
    return df.reset_index(drop=True)
