from __future__ import annotations
import copy
import re
from pathlib import Path
import yaml

DEFAULT_CONFIG_PATH = Path(__file__).resolve().parent.parent / "config.yaml"


def load_config(path: str | Path | None = None) -> dict:
    p = Path(path) if path else DEFAULT_CONFIG_PATH
    with open(p, "r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f)
    cfg = _coerce_numbers(cfg)
    cfg["_config_path"] = str(p)
    cfg["_root"] = str(p.parent)
    return cfg


_NUM = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")


def _coerce_numbers(obj):
    """PyYAML (YAML 1.1) reads '1.0e5' as a string; convert such strings to float recursively."""
    if isinstance(obj, dict):
        return {k: _coerce_numbers(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_coerce_numbers(v) for v in obj]
    if isinstance(obj, str) and _NUM.match(obj.strip()) and any(ch in obj for ch in "eE."):
        try:
            return float(obj)
        except ValueError:
            return obj
    return obj


def cfg_get(cfg: dict, dotted: str, default=None):
    cur = cfg
    for k in dotted.split("."):
        if not isinstance(cur, dict) or k not in cur:
            return default
        cur = cur[k]
    return cur


def deep_update(base: dict, upd: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in upd.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)
        else:
            out[k] = v
    return out
