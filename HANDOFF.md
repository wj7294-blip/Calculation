# 인수인계: 음극 SAM 연구 공백 조사 및 전지별 아이디어 (Robin 방식)

작성일: 2026-10-31. 이 문서는 이전 Claude Code 세션의 작업을 다른 세션에서 이어가기 위한 요약이다.
세션 링크: https://claude.ai/code/session_01SabbWXBmYwH1t2nGvwSLTi

## 1. 사용자가 요청한 것

"SAM(자기조립단분자막)이 음극에서도 아직 빈 부분이 있는지 봐줘. 리튬금속, 실리콘 음극, 아연전지, 다가이온전지 등 다각도로.
선행연구 조사를 많이 하고, 전지 시스템별로 아이디어를 내줘. 방식은 올해(2026) Nature에 나온 에이전트 기반 연구법을 참고."

- "올해 Nature 에이전트 연구법" = FutureHouse **Robin**, "A multi-agent system for automating scientific discovery",
  Nature 655, 497–505 (2026), DOI 10.1038/s41586-026-10652-y (arXiv 2505.13400, 코드 github.com/Future-House/robin).
- 사용자 배경: 이 저장소는 MACE uMLIP로 Cu(111)/Ag(111) 위 Li 흡착 에너지를 계산하는 프로젝트(`README.md`, `adsorption.py`, `verify.py`, `crosscheck.py`, `diagnose.py`).
  anode-free Li 금속 전지 맥락. 아이디어 제안 시 이 전산 파이프라인과 연결하도록 했다.

## 2. 산출물 위치 (모두 main 브랜치에 커밋·푸시됨, 최신 커밋 a2706ca)

| 파일 | 내용 |
|---|---|
| `reports/음극 SAM 연구 공백과 전지별 아이디어.md` | **최종 보고서(한국어, 약 97k자, 표 10개, 인용 약 170개)**. 사용자에게 전달 완료. |
| `research_notes/음극 SAM 연구 공백과 전지별 아이디어/robin_methodology.md` | Robin 워크플로 재구성(공개 코드 기반) + 부록 "Robin-by-hand recipe" |
| `.../sam_li_metal.md` | Li 금속 / anode-free Cu 위 SAM 선행연구 + 공백 표 |
| `.../sam_silicon.md` | Si 음극 SAM/grafting 선행연구(18행 부록 표) + 공백 |
| `.../sam_zinc.md` | Zn 음극 SAM(약 15–20편) + 포화/공백 지도 |
| `.../sam_multivalent_na_k.md` | Mg/Ca/Al/Na/K 금속 음극 SAM(사실상 n=1) + 널 결과 목록 |
| `.../sam_mechanisms_modeling.md` | 횡단 메커니즘(환원 탈착, 터널링, 쌍극자), 전산(DFT/MLIP), 품질관리 공백 |
| `.../sam_frontier_and_cathode_contrast.md` | 2025–26 최신 동향, 양극 측 SAM과 대조, 특허 |

보고서 구조: 1 Robin 손 적용 절차 → 2 시스템별 선행연구·공백 지도(2.1 Li, 2.2 Si, 2.3 Zn, 2.4 Mg/Ca/Al, 2.5 Na/K, 2.6 흑연·집전체·고체전해질, 2.7 횡단 요약표)
→ 3 횡단 공백 → 4 가설 후보 37개(CANDIDATE/HYPOTHESIS/REASONING, 시스템별 순위표; ID L1–L8, S1–S6, Z1–Z6, M1–M7, N1–N5, G1–G5)
→ 5 상위 아이디어 검증 실험·계산 설계(assay) → 6 통합 Top 5 → 7 한계·후속 확인 목록 → 결론.

## 3. 핵심 결론 요약

- 검증된 1차 논문 수: Li/anode-free Cu <10편, Si ~15편(대부분 공유결합 grafting), Zn 15–20편, Mg 1편(in situ 실란 SAM, Fan et al. Angew. 2025, 10.1002/anie.202416582), **Ca·Al·Na·K 금속 음극 및 Mg/Ca/Na/K용 집전체: 0편**.
- 횡단 공백 3가지: (a) 도금 전위(0 V vs Li+/Li)에서 SAM의 운명(잔존/탈착/무기물 전환)을 operando·전위 분해 XPS/ToF-SIMS로 추적한 논문 없음;
  (b) SAM–금속–알칼리금속 3성분 계면 MLIP/uMLIP 연구 0편; (c) "쌍극자 효과" 주장은 DFT 그림뿐, KPFM/UPS 측정 없음. 품질관리는 접촉각 단독.
- 포화 영역(피할 것): ZnSO4 전해질·두꺼운 Zn foil 위 단일 성분 소수성 SAM 대칭셀 수천 시간 보고.
- 통합 Top 5:
  1. L1+L2: 동일 Cu·전해질에서 head-out형 SAM(MPA, 5-methyltetrazole) vs 희생형 SAM(1,3-benzenedisulfonyl fluoride) 정면 비교 + 동위원소 표지·전위 분해 XPS/ToF-SIMS로 운명 추적
  2. L7(+Z4, G1): 반응성 없는 등구조(isostructural) 쌍극자 SAM 시리즈 + KPFM으로 "쌍극자 vs LiF 화학" 분리
  3. Z3: KOH Zn–air/Zn–Ni용 alkanethiol SAM(완전 공백, 인접 근거 강함)
  4. M1+M2: Mg 위 최초 ex situ phosphonic acid SAM, Cl-free vs Cl-rich 전해질 대비
  5. N1: anode-free Na용 Cu/Al 집전체 azole/MPA SAM (uMLIP에서 Li→Na adatom 치환으로 사전 선별)
- uMLIP 확장 제안(보고서 5절): thiolate/tetrazolate/phosphonate on Cu(111)/Ag(111) 1/3 ML를 PBE-D3와 비교 → SAM 덮인 Cu 위 Li adatom(기존 bare-Cu 기준 −0.42…−0.56 eV와 비교)
  → 환원 탈착 대리 지표 ΔE(Cu–S → Li–S–R) → NEB → Na/Mg adatom, Zn/Mg(0001) 슬랩. 3모델(MP-0b3/MPA-0/OMAT-0) 스프레드를 불확실성으로 사용.
  주의: 유기 단분자막–배터리 금속에 대한 uMLIP 검증 사례 없음, MACE-OFF는 금속 미포함.

## 4. 방법 (Robin 모방)

Robin 공개 코드에서 재구성한 절차: 5 쿼리 → 10 assay 아이디어 → 45 쌍대 비교 → Bradley–Terry(choix, α=0.1) → 최상위 assay;
10 쿼리 → 30 후보(CANDIDATE/HYPOTHESIS/REASONING) → 후보별 심층 보고서 → 300 무작위 쌍대 비교 → BT;
데이터 분석 5회 병렬 + 합의 → {tested, summary, questions, mechanistic insights} 4필드를 다음 라운드 프롬프트에 주입.
판정 기준 우선순위: 근거 강도 > 메커니즘 명확성 > 안전성 > 실현가능성 > 신규성. **순위는 Elo가 아니라 BT**.
이번 작업에서는 리서치 노트가 Crow/Falcon 역할, 순위는 LLM 토너먼트가 아닌 **수작업 H/M/L 적용**(보고서에 명시). 37개 후보 블록은 실제 LLM 쌍대 비교에 그대로 넣을 수 있는 형식.

## 5. 알려진 한계 (반드시 인지)

- 세션 네트워크 정책으로 Wiley·ACS·Elsevier·Nature·Science·RSC·PubMed·arXiv·Google Patents 등 전문 접근이 403으로 차단됨.
  **모든 수치는 검색 스니펫·초록 수준**. 노트와 보고서의 [UNVERIFIED], [PAPER-via-secondary], [memory, unverified] 표지를 그대로 유지했다.
- "0편" 진술은 노트당 30–45회 웹 검색 기준. Scopus/WoS/CNKI 전수 조사 아님. "molecular layer", "organic monolayer", "surface grafting" 등 다른 용어로 발표된 논문이 있을 수 있음.
- Robin 논문 원문·SI는 읽지 못했고 공개 코드 + 2차 요약으로 재구성.
- 리뷰들의 명시적 open-question 절은 읽지 못해 3절 횡단 공백은 1차 연구 커버리지로부터의 추론.

## 6. 다음 세션에서 할 만한 일 (우선순위 순)

1. **원문 검증**: 보고서 7절 목록. 특히 Gao et al. Nat. Energy 2020 (10.1038/s41560-020-0640-7)의 DFT 내용, Liang et al. Adv. Mater. 2026 (10.1002/adma.73830)의 600 vs 505 Wh/kg 불일치,
   Electrochim. Acta 2024 MPA/MUA·HMDS 논문의 전해질·풀셀 데이터, Fan 2025의 CE·전해질, Science 2022 abn1818 저자·수치, Robin Nature 원문 Methods/SI.
   → 전문 접근이 되는 환경(기관 VPN, 또는 네트워크 허용 도메인 추가)에서 수행.
2. **실제 Robin 토너먼트 실행**: 보고서 4절의 37개 CANDIDATE 블록을 LLM 판정자(Robin 기준 순서 프롬프트)로 300 쌍대 비교 후 `choix.ilsr_pairwise`(α=0.1)로 BT 점수화. 수작업 순위와 비교.
3. **uMLIP 사다리 1단**: 저장소 파이프라인(`adsorption.py` 구조 재사용)으로 methanethiolate / tetrazolate / methylphosphonate on Cu(111) 1/3 ML 흡착 에너지·기하를 MP-0b3/MPA-0/OMAT-0로 계산, 가능하면 PBE-D3 단일점과 비교.
   선행 과제: `results/diagnosis/structures/` 구조로 PBE 단일점을 수행해 Ag–Cu Li 결합 부호 문제(README "원인 점검" 절)를 먼저 닫을 것.
4. Scopus/WoS에서 Na/K/Ca/Al + "monolayer|molecular layer|grafting|anode" 재검색으로 "0편" 진술 재확인.
5. 특허: Samsung SDI·LGES·SK On·CATL·QuantumScape·SES의 SAM + anode/current collector 청구항, US 12,658,447 양수인 확인.

## 7. 작업 관례

- 커밋 메시지 끝에 `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` 및 세션 링크를 붙였음(세션 지침). 새 세션은 그 세션의 지침을 따를 것.
- 저장소 stop hook이 미커밋 파일이 있으면 커밋·푸시를 요구함(`~/.claude/stop-hook-git-check.sh`). main 브랜치에 직접 커밋해 왔음.
- 리서치 노트 형식: Key Question별 Takeaway / Cited Findings(인라인 URL) / Inferences / Gaps.
