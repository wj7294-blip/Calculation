# Calculation

adsorption energy calculation

사전학습 범용 머신러닝 포텐셜(uMLIP) MACE로 Cu(111) 표면 위 Li 원자 1개의 흡착 에너지를 자리별로 계산한다. GPU 없이 CPU만으로 실행한다.

## 실행 방법

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install torch --index-url https://download.pytorch.org/whl/cpu   # CPU 전용 PyTorch
pip install -r requirements.txt
python adsorption.py
```

- 처음 실행할 때 MACE 체크포인트(약 76 MB)를 GitHub release에서 `~/.cache/mace/`로 자동으로 받는다.
- 표면 금속, 흡착 원자, 슬랩 크기, 층수, 진공 두께, 고정 층수, 수렴 기준, 모델은 `adsorption.py` 맨 위 변수에서 바꾼다.

> 참고: 이번 실행 환경(Claude Code 클라우드 세션)에서는 `download.pytorch.org`가 네트워크 정책으로 막혀 있었다. 그래서 PyPI 기본 `torch` 휠(`2.14.1+cu130`)을 설치했다. GPU가 없으니 CUDA는 쓰지 않았고(`torch.cuda.is_available() == False`) 계산은 모두 CPU에서 했다.

## 계산 조건

| 항목 | 값 |
|---|---|
| 모델 | MACE-MP-0b3 medium (`mace_mp(model="medium-0b3")`, [mace_mp_0b3 release](https://github.com/ACEsuit/mace-mp/releases/tag/mace_mp_0b3)) |
| 소프트웨어 | mace-torch 0.3.16, torch 2.14.1, ase 3.29.0, e3nn 0.4.4, float64 |
| 하드웨어 | CPU 4 스레드, GPU 없음 |
| 표면 | Cu(111), 4×4, 4층(64 원자), 진공 15 Å(ASE `fcc111` 기준 위아래 각각), 아래 2층 고정 |
| 격자상수 | MACE로 최적화한 Cu fcc a = 3.614 Å |
| 흡착 원자 | Li 1개, 초기 높이 2.0 Å |
| 흡착 자리 | fcc hollow, hcp hollow, top, bridge |
| 최적화 | BFGS, fmax = 0.05 eV/Å (bulk는 셀까지 포함해 0.01 eV/Å) |
| 기준 에너지 | 깨끗한 슬랩(같은 조건으로 최적화), Li bcc bulk 원자당 에너지(a = 3.417 Å, −1.9008 eV/atom) |

E_ads = E(slab+Li) − E(slab) − E(Li bulk, 원자당)

E_ads가 음수이면, Li가 bulk 금속 Li로 있을 때보다 Cu 표면에 흡착했을 때 더 안정하다는 뜻이다.

## 결과

| 흡착 자리 | E_ads (eV) | 표면에서 Li까지 높이 (Å) | BFGS steps | 계산 시간 (s) |
|---|---:|---:|---:|---:|
| fcc hollow | −0.537 | 2.061 | 4 | 5.5 |
| **hcp hollow** | **−0.563** | 2.051 | 4 | 4.1 |
| top | −0.364 | 2.022 | 12 | 11.1 |
| bridge | −0.527 | 2.068 | 5 | 4.8 |

- 높이는 최적화 후 Li의 z 좌표에서 맨 위 Cu 층 평균 z를 뺀 값이다.
- 자리 간 안정성 순서는 hcp ≳ fcc > bridge ≫ top이다. hcp와 fcc의 차이(약 26 meV)는 모델 오차 범위 안일 수 있다.
- 최적화하는 동안 Li의 수평 이동은 모두 0.01 Å 이하였다. top과 bridge는 대칭 때문에 수평 힘이 0이라 그 자리에 머문 것이다. 실제로 안정한 극소점이 아니라 hollow 자리 사이의 안장점일 수 있다(진동수 분석으로는 확인하지 않았다).

**계산 시간**(CPU 4 스레드): bulk 기준 1.3 s, 깨끗한 슬랩 1.5 s(초기 힘이 이미 fmax보다 작아 0 step), 흡착 4개 자리 합계 25.5 s다. 모델 로딩까지 포함한 전체 실행은 약 38 s다(모델 다운로드 시간 제외).

![E_ads by site](results/adsorption_energies.png)

### 결과 파일 (`results/`)

| 파일 | 내용 |
|---|---|
| `adsorption_energies.csv` | 자리별 E_ads, 높이, 수평 이동, step 수, 시간, 총에너지 |
| `adsorption_energies.png` | 자리별 E_ads 막대그래프 |
| `Cu111_Li_<site>.traj`, `clean_slab.traj` | BFGS 최적화 궤적(ASE trajectory, 마지막 프레임이 최종 구조) |
| `POSCAR/POSCAR_<name>` | 최종 구조(VASP POSCAR, direct 좌표) |
| `logs/*.log` | BFGS 로그 |
| `run_info.json` | 모델, 패키지 버전, 설정, 기준 에너지 |

## 한계

- **사전학습 uMLIP 결과는 경향 비교용이다.** MACE-MP-0b3는 Materials Project(PBE/PBE+U) bulk 데이터 위주로 학습했다. 표면이나 흡착 에너지의 절대값에는 수십~수백 meV 수준의 오차가 있을 수 있다. 정확한 값이 필요하면 같은 구조(`results/POSCAR/`)로 DFT 계산을 해서 검증해야 한다.
- **단일 Li 원자 흡착이다.** 4×4 셀에 Li 1개(커버리지 1/16 ML)만 넣었기 때문에, Li 사이 상호작용이나 커버리지가 높아질 때의 효과는 반영되지 않았다.
- 슬랩 4층에 아래 2층을 고정했고, 영점 에너지와 온도 효과는 넣지 않았다. top과 bridge 자리의 안정성(극소점인지 안장점인지)도 검증하지 않았다.
