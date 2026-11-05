# echem_tools

전기화학 실험 데이터 분석 도구 (Python, CPU 전용).

| 모듈 | 기능 |
|---|---|
| **EIS → DRT** | 전처리(주파수 자르기, 유도성 꼬리), lin-KK 검증, DRT(가우시안 RBF + Tikhonov, 비음수, GCV λ), 피크 표, Nyquist/Bode/DRT 그림, 비교 모드 |
| **고전압 플로팅** | 정전압 구간 자동 검출, 누설전류(mA/g, µA/cm²), 누적 누설 용량(전체 / 과도구간 제외), log-log 감쇠 기울기, 그룹 비교(이원 ANOVA + 사후검정 + 비모수), 상호작용 그림 |
| **공통** | BioLogic `.mpt`, Gamry `.DTA`, CSV(일반/Neware/LANHE, `config.yaml`로 열 매핑) 로더, `metadata.csv` 기반 셀 관리, `results/` 에 CSV·PNG·`report.html` |
| **검증** | 등가회로 합성 스펙트럼 + 합성 플로팅 데이터(5 표면 × 3 전해질 × n=4, 상호작용 심어둠)로 전체 파이프라인 자동 검증 |

---

## 1. 설치

```bash
git clone https://github.com/wj7294-blip/Calculation.git
cd Calculation/echem_tools
pip install -r requirements.txt
```

Python ≥ 3.10. `pyDRTtools` 는 선택 사항입니다(설치돼 있고 import 가 되면 `eis.drt.backend: auto` 에서 먼저 시도하고, 실패하면 자동으로 내장 솔버를 씁니다. 현재 배포판은 PyQt5 를 강제로 import 하므로 headless 환경에서는 대부분 내장 솔버가 사용됩니다).

## 2. 빠른 시작 (파일만 올리면 돌아감)

### 2-1. EIS → DRT: 파일을 직접 지정

```bash
python analyze.py drt cell01.mpt cell02.DTA cell03.csv
```

```bash
python analyze.py drt data/eis/*.mpt --lambda 1e-3 --inductive drop --fmin 0.1 --out results/eis_run1
```

- 파일을 여러 개 주면 각각 분석하고 **비교 그림**(`results/eis/compare_*.png`)도 만듭니다.
- `--meta metadata.csv` 를 주면 파일 이름을 metadata 의 `file` 열과 맞춰 `cell_id` 라벨과 그룹(`--group-by surface` 등)을 붙입니다.
- `--lambda auto`(기본, GCV) 또는 숫자. `--inductive fit|drop|keep`.
- 파일 없이 `python analyze.py drt --meta metadata.csv` 만 주면 metadata 의 `type=eis` 행을 모두 분석합니다.

### 2-2. 플로팅 누설전류: metadata 기반 그룹 분석

```bash
python analyze.py floating --meta metadata.csv
```

```bash
python analyze.py floating --meta metadata.csv --data-dir data --transient-h 2 --out results/float_run1
```

단일 셀을 빨리 보고 싶을 때(그룹 통계 없음):

```bash
python analyze.py floating cell07.csv --mass-mg 8.1 --area-cm2 1.13
```

### 2-3. 한 줄 전체 실행 / 합성 데이터 검증

```bash
python analyze.py all --meta metadata.csv          # metadata 의 eis + float 행 전부 → results/report.html
```

```bash
python analyze.py validate                         # 합성 데이터 생성(data_synthetic/) + 전체 검증(results_synthetic/report.html)
```

### 2-4. Colab

`colab_analysis.ipynb` 를 Colab 에서 열고 순서대로 실행: ① 설치(기본값: `wj7294-blip/Calculation` 저장소의 `echem_tools/` 폴더를 clone) → ② 파일 업로드(여러 개, `metadata.csv` 포함 가능) → ③ 분석 종류(DRT / Floating)와 설정값(λ, 과도구간 시간, V_hold, 질량/면적, CSV 프로파일) → ④ 실행하면 표·그림이 노트북에 표시되고 `report.html` 과 결과 zip 이 다운로드됩니다.

## 3. 파일 넣는 위치와 metadata 작성법

```
echem_tools/
├─ data/
│   ├─ eis/     ← EIS 파일 (.mpt / .DTA / .csv)
│   └─ float/   ← 플로팅(CC-CV) 파일
├─ metadata.csv
├─ config.yaml
└─ results/     ← 출력
```

`metadata.csv` (쉼표 구분, `#` 로 시작하는 줄은 무시):

```
cell_id,file,type,surface,electrolyte,active_mass_mg,area_cm2,note
B-A-1,float/B-A-1.mpt,float,bare,A,8.02,1.13,
S8-TTE-1,float/S8-TTE-1.csv,float,SAM8,TTE,7.95,1.13,Neware
E-01,eis/E-01.DTA,eis,bare,A,,,"after 10 cycles"
```

| 열 | 의미 |
|---|---|
| `cell_id` | 고유 라벨(그림·표에 표시). 중복 불가 |
| `file` | `data/` 기준 상대경로(또는 절대경로) |
| `type` | `eis` 또는 `float`. 비우면 폴더 이름(eis/float)으로 추정 |
| `surface`, `electrolyte` | 그룹 비교 요인(플로팅). EIS 에서는 비교 모드 그룹 키로 사용 가능 |
| `active_mass_mg`, `area_cm2` | 정규화에 사용(없으면 해당 지표 NaN) |
| `note`(또는 `비고`) | 자유 메모 |
| `profile` (선택) | CSV 열 매핑 프로파일(`generic`/`neware`/`lanhe`). 없으면 자동 시도 |

### CSV 열 이름 매핑 (`config.yaml` → `csv_profiles`)

Neware/LANHE 내보내기 열 이름이 다르면 `candidates` 목록에 추가하면 됩니다. `scale` 은 단위 변환(예: A→mA = 1000).

```yaml
csv_profiles:
  neware:
    time_s:     {candidates: ["Total Time(s)", "TestTime(s)"], scale: 1.0}
    voltage_V:  {candidates: ["Voltage(V)"], scale: 1.0}
    current_A:  {candidates: ["Current(A)"], scale: 1000.0}   # A -> mA
    step:       {candidates: ["Step Type"], scale: 1.0}
    cv_step_labels: ["CV_Chg", "CCCV_Chg"]                     # 정전압 구간 판별 보조
```

프로파일을 지정하지 않은 CSV 는 `default_csv_profile` → 나머지 프로파일 순으로 자동 시도합니다.

### 지원 형식 세부

- **BioLogic .mpt**: `Nb header lines` 로 헤더 분리, 소수점 쉼표 처리. EIS: `freq/Hz, Re(Z)/Ohm, -Im(Z)/Ohm`; 크로노: `time/s, Ewe/V, <I>/mA (또는 I/mA)`.
- **Gamry .DTA**: `ZCURVE TABLE`(EIS: `Freq, Zreal, Zimag`), `CURVE*/OCVCURVE TABLE`(크로노: `T, Vf, Im`[A→mA 변환]).
- 전류 부호: **충전(산화) = 양수**. 시간이 스텝마다 0 으로 리셋되는 내보내기는 자동으로 누적 시간으로 고칩니다.

## 4. 출력물

`results/` (또는 `--out`):

| 파일 | 내용 |
|---|---|
| `report.html` | 모든 표·그림(이미지 내장) 요약 리포트 |
| `eis_summary.csv` | 스펙트럼별 R∞, L, R_pol, λ, 적합 잔차, lin-KK 잔차/경고 |
| `eis_peaks.csv` | 피크별 τ, f=1/(2πτ), 저항(면적), R_pol 대비 비율, 적분 구간 |
| `eis/<label>_drt.csv / _fit.csv / _drt.png / _gcv.png / _kk.png` | γ(τ), 적합 Z, Nyquist·Bode·DRT 그림, GCV/L-curve, KK 잔차 |
| `eis/compare_<group>.png` | 비교 모드 겹쳐 그리기 |
| `float_cells.csv` | 셀별 지표(아래 정의) + 이상치 플래그 |
| `float_group_summary.csv` | 그룹별 평균/SD/n |
| `float_anova.csv`, `float_posthoc_tukey.csv`, `float_nonparametric.csv` | 이원 ANOVA, Tukey, SRH/Kruskal/Dunn |
| `float/<cell>_cv.csv / .png`, `float/all_cells_decay.png`, `float/group_*.png`, `float/interaction_*.png` | 셀별 정전압 구간 데이터·그림, 그룹 비교·상호작용 그림 |

---

## 5. 지표 정의와 계산식

### 5-1. EIS 전처리 / lin-KK

- 주파수 범위 `[freq_min_Hz, freq_max_Hz]` 밖은 제거. 중복 주파수 제거.
- 유도성 꼬리(−Z'' < 0 인 고주파 점): `inductive_tail: fit`(기본; 점은 유지하고 인덕턴스 L 을 함께 적합) / `drop`(점 제거) / `keep`.
- **lin-KK** (Schönleber 2014, `impedance.py` 의 구성 함수 사용): τ 를 측정 주파수 범위에 로그 등간격으로 고정한 M 개 RC + R + L + C 를 선형 최소제곱으로 적합. M 은 3 부터 늘리며 μ = 1 − Σ|R_k<0| / Σ|R_k≥0| < c(0.85) 를 처음 만족하는 M 을 택하되, `impedance.linKK` 가 너무 작은 M 에서 멈추는 문제를 막기 위해 **rms 잔차가 최소값의 `plateau_tol`(1.2)배 이내**인 가장 작은 M 을 고릅니다.
  잔차 Δ = (Z_meas − Z_KK)/|Z_meas| 의 **rms > 1 %** 또는 **max > 5 %** 이면 경고(`kk_warn=True`). 경고는 드리프트·비선형·노이즈 가능성을 뜻하며 DRT 해석에 주의가 필요합니다.

### 5-2. DRT

모델 (γ 는 ln τ 에 대한 분포):

```
Z(ω) = R∞ + jωL + ∫ γ(ln τ) / (1 + jωτ) d ln τ ,   γ(ln τ) = Σ_m x_m φ_m(ln τ)
φ_m(y) = exp(−(ε (y − ln τ_m))²)   (가우시안 RBF, τ 격자 10/decade, 측정 범위 ±1 decade 확장)
```

- 기저 폭: FWHM = `rbf_shape`(1.5) × 격자 간격(ln 단위). 1 보다 작으면 기저가 격자보다 좁아 γ 가 울퉁불퉁해집니다.
- 풀이: `min ‖W(Ax − b)‖² + λ‖Ux‖²  s.t. x ≥ 0` (NNLS). W = 1/|Z|(`weighting: modulus`) 또는 1. UᵀU = ∫ φ⁽ᵏ⁾ᵢ φ⁽ᵏ⁾ⱼ d ln τ (k = `reg_order`, 기본 2차 미분 Tikhonov). **R∞ ≥ 0, L ≥ 0 은 벌점 없이 같은 문제의 별도 열**로 적합됩니다.
- λ 는 **무차원**으로 정규화했습니다: tr(UᵀU) = tr(AᵀWᵀWA) 가 되도록 U 를 스케일 → λ ≈ 1 이면 데이터 적합과 매끄러움을 같은 비중으로 봄. 자동 선택은 비제약 Tikhonov 해에 대한 GCV
  `GCV(λ) = n‖Ax_λ − b‖² / (n − tr H_λ)²` 의 최소값(`lambda_grid` 1e-6 ~ 10, 36 점). `lambda: 1e-3` 처럼 숫자를 주면 수동.
- 피크: γ(τ)에서 `scipy.signal.find_peaks`(prominence ≥ 3 % max γ). 피크 경계는 이웃 피크 사이 최소점, **피크 저항 R_i = ∫_경계 γ d ln τ**, f_i = 1/(2πτ_i). R_pol = Σ x_m √π/ε = ∫γ d ln τ.

### 5-3. 플로팅 누설전류

프로토콜: 정전류(C-rate `c_rate`)로 `V_hold`(4.6 V) 까지 충전 → `V_hold` 정전압 `hold_hours`(24 h) 유지, I–t 기록.

- **정전압 구간 검출**: 전압이 처음 `V_hold − V_tol` 에 도달한 시점부터, `|V − V_hold| ≤ V_tol`(0.02 V) 이고 I > 0 인 연속 구간(3점 초과 이탈 시 종료). Neware/LANHE CSV 에 `cv_step_labels` 에 맞는 step 라벨이 있으면 그 구간(가장 긴 연속 구간)을 우선 사용. 검출 방식은 `cv_detect` 열에 기록.
- **누설전류** `leak_mA_g`, `leak_uA_cm2`: 정전압 구간 **마지막 `leak_window_h`(1 h) 의 평균 전류** / 활물질 질량(g), × 1000 / 면적(cm²).
- **누적 누설 용량** `Q_leak_mAh_g` = ∫₀ᵀ I dt / 3600 / m (정전압 구간 전체, 사다리꼴 적분); `Q_leak_ex_transient_mAh_g` = ∫_{t_tr}ᵀ I dt / 3600 / m, t_tr = `transient_h`(1 h, 설정 가능). 차이 `Q_transient_mAh_g` 는 초기 과도(이중층 충전 + 농도 완화 + 초기 분해)에 해당.
- **감쇠 기울기** `decay_slope`: `decay_fit_window_h`(2–24 h) 구간에서 log₁₀ I 를 log₁₀ t 에 선형 회귀한 기울기(−n). 순수 확산/과도이면 ≈ −0.5, 정상상태 누설이 지배하면 0 에 가까워짐(기울기의 절대값이 작을수록 누설 전류가 꾸준히 흐름). `decay_r2`, 초기 1 h 기울기 `decay_slope_early` 도 기록.
- **그룹 통계** (`surface × electrolyte`, 그룹당 n≈4):
  - 그룹별 평균 ± SD, 개별 셀 점, 이상치는 그룹 내 1.5·IQR(`outlier_method: iqr`, 또는 MAD) 로 **표시만**(제외하지 않음, `outlier_<metric>` 열).
  - 이원 분산분석(type II SS, 상호작용 포함, `statsmodels`), 부분 η². 가정 검정: 잔차 Shapiro–Wilk, 그룹 Levene(중앙값).
  - 사후검정: Tukey HSD — (i) 전체 (surface|electrolyte) 조합 쌍, (ii) **각 전해질 내 surface 쌍**, (iii) 각 요인 주효과.
  - 정규성/등분산이 의심되면(`nonparametric: auto`; `always` 로 강제 가능) 순위 기반 **Scheirer–Ray–Hare**(이원, 상호작용 포함), Kruskal–Wallis(요인별), Dunn(Holm) 사후검정을 함께 보고.
  - 상호작용 그림: x = surface, 선 = electrolyte. 선이 평행하지 않으면(상호작용 p < 0.05) "SAM 효과가 전해질에 따라 다르다"는 뜻입니다. 순위 기반 검정은 순서가 바뀌지 않는(비교차) 상호작용에는 둔감하다는 점에 유의.

## 6. DRT λ 선택이 결과에 미치는 영향

`results_synthetic/eis/lambda_sweep.png` (τ 비 8 인 두 RC, 1 % 노이즈):

- **λ 가 작으면**(≲1e-5): 피크가 날카롭고 저항 적분은 정확하지만, 노이즈가 **가짜 작은 피크**나 **하나의 피크가 둘로 갈라지는** 형태로 나타날 수 있음. 피크 수를 그대로 믿으면 안 됨.
- **λ 가 크면**(≳1e-1): 피크가 넓어지고 낮아지며, 가까운 피크(τ 비 ≲ 10)가 **하나로 합쳐짐**. 합쳐진 피크의 면적(저항 합)은 대체로 보존되지만 τ 위치는 가중 평균으로 이동.
- GCV 는 보통 1e-4 ~ 1e-3 부근을 고르며, 합성 검증에서 피크 τ 오차 < 0.03 decade, RC 저항 오차 수 % 수준이었습니다. 단, GCV 는 노이즈가 KK 를 위반할 때(드리프트) 지나치게 작은 λ 를 고를 수 있으므로, **같은 시료 시리즈를 비교할 때는 λ 를 하나로 고정**(`--lambda 1e-3`)하는 편이 안전합니다.
- 피크 분리 한계: 기저 폭과 λ 때문에 τ 비 3 이하인 과정은 사실상 분리되지 않습니다. 또한 유한길이 Warburg 같은 분포형 요소는 DRT 에서 단일 피크가 아니라 τ_k = 4τ_W/((2k−1)²π²), R_k = 8R_W/((2k−1)²π²) 급수(k=1 항이 81 %)로 나타나므로, 저주파 주 피크 옆의 작은 피크들이 인접 RC 피크와 합쳐져 그 RC 저항을 과대평가할 수 있습니다(검증 표의 RC3 +14 % 가 그 예).

## 7. 한계

- DRT 는 τ 격자 범위 밖(측정 주파수 ±1 decade)의 과정을 표현하지 못하며, 반무한 Warburg(ω^-1/2)·CPE 처럼 ln τ 에서 끝이 없는 분포는 γ 가 격자 끝에서 잘려 R_pol 이 과소평가됩니다. 저주파 끝에서 γ 가 올라가 있으면 측정을 더 낮은 주파수까지 확장하세요.
- lin-KK 경고가 있는 스펙트럼의 DRT 는 정량 해석보다 정성 비교에만 쓰는 것이 안전합니다.
- 누설전류 "마지막 1 h 평균"에는 24 h 시점에도 남아 있는 과도 전류 꼬리(t^-n)가 포함되어 정상상태 누설보다 체계적으로 큽니다(합성 검증에서 +10 ~ +15 %). 그룹 간 비교에는 영향이 작지만 절대값으로 해석할 때 유의하세요. 더 긴 유지 시간이나 두 지수 모델 적합이 필요하면 `decay_fit_window_h` 를 조정하거나 `float/<cell>_cv.csv` 를 직접 쓰세요.
- 정전압 구간 검출은 전압 창 기준이므로 IR 강하가 큰 셀에서 CC→CV 전이 수 점이 포함되거나 빠질 수 있습니다(`V_tol` 조정). 전류 부호가 반대(충전 음수)인 내보내기는 로더 또는 `scale: -1` 로 뒤집어야 합니다.
- 그룹당 n=4 는 정규성 검정의 검정력이 매우 낮습니다. `nonparametric: always` 로 두 결과를 항상 함께 보는 것을 권장합니다. 순위 기반 이원 검정(SRH)은 근사 검정입니다.
- pyDRTtools 백엔드는 설치 환경에 따라 import 가 되지 않을 수 있고(PyQt5 의존), 그 경우 자동으로 내장 솔버로 대체됩니다. 두 솔버의 λ 는 서로 다른 단위이므로 직접 비교하지 마세요.

## 8. 검증 결과 요약 (`python analyze.py validate`)

합성 회로 (R0 + 2RC + L), (R0 + 3RC + 유한길이 Warburg + L), (R0 + 가까운 2RC), 노이즈 0.5 / 1 %, 세 로더(.mpt/.DTA/.csv) 왕복:

| 항목 | 결과 |
|---|---|
| RC 피크 τ 오차 | < 0.03 decade |
| RC 피크 저항 오차 | 0.1–3.5 % (Warburg 급수와 겹치는 RC3 는 +14–16 %) |
| R∞ 오차 | < 1 % |
| FLW 주 피크 (k=1) | τ 오차 < 0.02 decade, 저항 오차 ≈ +4 % |
| 플로팅: 셀 누설전류 복원(노이즈 없는 모델의 마지막 1 h 평균 대비) | 평균 0.1 %, SD 0.9 % |
| 이원 ANOVA | surface, electrolyte, 상호작용 모두 p < 1e-4 (심어둔 SAM8|TTE 시너지, SAM10|GSE·SAM9|GSE 길항 검출) |

전체 표와 그림은 `results_synthetic/report.html` 참조.

## 9. 테스트

```bash
python -m pytest -q tests
```

## 10. 구조

```
analyze.py                  CLI (drt / floating / all / validate)
colab_analysis.ipynb        Colab 노트북
config.yaml                 모든 설정
echem/
  config.py                 YAML 로드
  io/loaders.py             .mpt / .DTA / CSV 로더, 열 매핑
  io/metadata.py            metadata.csv
  eis/preprocess.py         자르기, 유도성 꼬리, lin-KK
  eis/drt.py                DRT 솔버(내장) + pyDRTtools 훅, 피크 검출
  eis/synth.py              등가회로 합성 스펙트럼
  eis/plots.py              Nyquist/Bode/DRT, GCV, λ sweep
  floating/analysis.py      정전압 구간, 누설전류, 누적 용량, 감쇠 기울기
  floating/stats.py         ANOVA, Tukey, SRH, Kruskal, Dunn, 이상치
  floating/synth.py         합성 플로팅 데이터 (.mpt/.DTA/Neware CSV)
  floating/plots.py         셀 곡선, 그룹 막대, 상호작용
  pipeline_eis.py, pipeline_float.py, report.py
tests/                      pytest
results_synthetic/          합성 검증 결과(report.html 포함)
```
