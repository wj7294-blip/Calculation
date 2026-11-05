"""전체 라이브러리 실행: library.csv → results/summary.csv, results/molecules/*.json, xyz/, molden/.

사용법:  python run_pipeline.py [--nconf 20] [--workers 4] [--only NAME ...]
         python run_pipeline.py --convergence   # 5개 분자 20 vs 50 conformer 검증
"""
import argparse
import csv
import json
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")

CONVERGENCE_SET = ["SH-Ph-CF3", "PO3H2-C3-OCH3", "COOH-C3-SO3H", "ZW-SO3-n3-m4", "L-cysteine-zwitterion"]


def load_library(path="library.csv"):
    rows = []
    for r in csv.DictReader(open(path)):
        for k in ("n", "m"):
            r[k] = int(r[k]) if r[k] else ""
        r["exp_ip"] = float(r["exp_ip"]) if r["exp_ip"] else ""
        rows.append(r)
    return rows


def _worker(args):
    row, n_conf, seed, outdir, tag = args
    from compute import process_molecule
    try:
        return process_molecule(row, n_conf=n_conf, seed=seed, outdir=outdir, tag=tag,
                                xyzdir="xyz" if not tag else "results/convergence/xyz",
                                moldendir="molden" if not tag else "results/convergence/molden",
                                do_ipea=not tag)
    except Exception as e:  # noqa
        import traceback
        return dict(name=row["name"] + tag, error=traceback.format_exc())


def rebuild_summary(outdir, summary_path):
    """outdir 의 모든 분자 JSON 에서 summary.csv 를 다시 만든다 (부분 재실행 후에도 전체 표 유지)."""
    import glob
    rows = [json.load(open(p))["summary"] for p in sorted(glob.glob(os.path.join(outdir, "*.json")))]
    if not rows:
        return
    cols = list(rows[0].keys())
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for s in sorted(rows, key=lambda s: s["name"]):
            w.writerow(s)


def log_line(msg, logf):
    print(msg, flush=True)
    logf.write(msg + "\n")
    logf.flush()


def run(rows, n_conf, workers, outdir, summary_path, tag="", seed=42):
    os.makedirs("logs", exist_ok=True)
    logf = open("logs/progress.log", "a")
    log_line(f"# {time.strftime('%Y-%m-%d %H:%M:%S')} 시작: 분자 {len(rows)}개, conformer {n_conf}개, workers {workers}", logf)
    jobs = [(r, n_conf, seed, outdir, tag) for r in rows]
    results, t0 = [], time.time()
    with Pool(workers) as pool:
        for i, s in enumerate(pool.imap_unordered(_worker, jobs), 1):
            if "error" in s:
                log_line(f"[{i}/{len(rows)}] {s['name']} | 실패: {s['error'].strip().splitlines()[-1]}", logf)
            else:
                log_line(f"[{i}/{len(rows)}] {s['name']} | conformer {s['n_conf_ok']}개 | 최저 E {s['E_min_Eh']:.6f} Eh"
                         f" | HOMO {s['HOMO']:.2f} eV | {s['time_s']:.1f} s", logf)
            results.append(s)
    log_line(f"# 완료: 총 {time.time()-t0:.0f} s", logf)
    rebuild_summary(outdir, summary_path)
    errs = [s for s in results if "error" in s]
    if errs:
        with open(summary_path.replace(".csv", "_errors.json"), "w") as f:
            json.dump(errs, f, indent=1)
    return results


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nconf", type=int, default=20)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--convergence", action="store_true")
    a = ap.parse_args()
    rows = [r for r in load_library() if not r["duplicate_of"]]
    if a.convergence:
        rows = [r for r in rows if r["name"] in CONVERGENCE_SET]
        os.makedirs("results/convergence", exist_ok=True)
        run(rows, 50, a.workers, "results/convergence", "results/convergence/summary_50.csv", tag="_50conf", seed=7)
    else:
        if a.only:
            rows = [r for r in rows if r["name"] in a.only]
        run(rows, a.nconf, a.workers, "results/molecules", "results/summary.csv")
