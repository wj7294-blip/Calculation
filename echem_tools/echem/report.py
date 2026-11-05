"""Minimal self-contained HTML report (images embedded as base64)."""
from __future__ import annotations
import base64
import datetime as _dt
from pathlib import Path
import pandas as pd
from jinja2 import Template

_TPL = Template("""<!doctype html><html lang="ko"><head><meta charset="utf-8">
<title>{{ title }}</title>
<style>
body{font-family:-apple-system,"Segoe UI",Roboto,"Noto Sans KR",sans-serif;margin:0;background:#f6f6f4;color:#0b0b0b}
main{max-width:1200px;margin:0 auto;padding:24px 16px}
h1{font-size:22px;margin:0 0 4px} h2{font-size:17px;margin:32px 0 8px;border-bottom:1px solid #d9d8d2;padding-bottom:4px}
h3{font-size:14px;margin:18px 0 6px;color:#52514e}
.meta{color:#52514e;font-size:12px;margin-bottom:16px}
table{border-collapse:collapse;font-size:12px;margin:8px 0;background:#fff} th,td{border:1px solid #e0dfda;padding:3px 7px;text-align:right}
th{background:#eeede8;text-align:center} td:first-child,th:first-child{text-align:left}
img{max-width:100%;height:auto;background:#fff;border:1px solid #e0dfda;margin:6px 0}
.warn{color:#b3261e;font-weight:600} .ok{color:#0a7a0a} .note{font-size:12px;color:#52514e}
pre{background:#fff;border:1px solid #e0dfda;padding:8px;font-size:11px;overflow:auto}
.scroll{overflow-x:auto}
</style></head><body><main>
<h1>{{ title }}</h1>
<div class="meta">generated {{ now }} · echem_tools {{ version }}</div>
{% for s in sections %}
<h2>{{ s.title }}</h2>
{% for b in s.blocks %}
  {% if b.kind == 'html' %}{{ b.html|safe }}
  {% elif b.kind == 'text' %}<p class="note">{{ b.text }}</p>
  {% elif b.kind == 'table' %}{% if b.title %}<h3>{{ b.title }}</h3>{% endif %}<div class="scroll">{{ b.html|safe }}</div>
  {% elif b.kind == 'image' %}{% if b.title %}<h3>{{ b.title }}</h3>{% endif %}<img src="data:image/png;base64,{{ b.b64 }}" alt="{{ b.title }}">
  {% elif b.kind == 'pre' %}<pre>{{ b.text }}</pre>{% endif %}
{% endfor %}
{% endfor %}
</main></body></html>""")


class Report:
    def __init__(self, title: str):
        self.title = title
        self.sections: list[dict] = []

    def section(self, title: str) -> "Report":
        self.sections.append({"title": title, "blocks": []})
        return self

    def _add(self, b):
        if not self.sections:
            self.section("Results")
        self.sections[-1]["blocks"].append(b)

    def text(self, t: str):
        self._add({"kind": "text", "text": t})

    def html(self, h: str):
        self._add({"kind": "html", "html": h})

    def pre(self, t: str):
        self._add({"kind": "pre", "text": t})

    def table(self, df: pd.DataFrame, title: str = "", fmt: str = "{:.4g}", max_rows: int = 200):
        if df is None or len(df) == 0:
            self._add({"kind": "text", "text": f"{title}: (empty)"})
            return
        d = df.head(max_rows).copy()
        html = d.to_html(index=not isinstance(d.index, pd.RangeIndex), float_format=lambda v: fmt.format(v), border=0, na_rep="")
        self._add({"kind": "table", "title": title, "html": html})

    def image(self, path: Path | str | None, title: str = ""):
        if path is None or not Path(path).exists():
            return
        b64 = base64.b64encode(Path(path).read_bytes()).decode()
        self._add({"kind": "image", "title": title, "b64": b64})

    def save(self, out: Path | str) -> Path:
        from . import __version__
        out = Path(out); out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(_TPL.render(title=self.title, now=_dt.datetime.now().strftime("%Y-%m-%d %H:%M"),
                                   version=__version__, sections=self.sections), encoding="utf-8")
        return out
