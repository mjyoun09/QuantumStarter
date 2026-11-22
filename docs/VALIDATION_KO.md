# DFT 검증 기록 · 최근 정리 2026-11-22

0.10.0에서는 UV–Vis 피크 24개를 PySCF와 대조하고, 세 기체의 열물성 108행과 NASA7 계수를 검증했다. 결합 해리 스캔은 H₂와 물에서 실수 스핀 곡률의 부호 전환과 더 낮은 UKS 에너지를 재현했다. [조건·수치·재현 명령](MOLECULAR_SUITE_KO.md)

## H₂ 결합 해리에서 RKS와 UKS가 갈라지는 거리

QS-20261121-H2-DISS. 자체 분자 PBE/def2-SVP로 H–H 거리 0.55–4.00 Å의 19점을 계산했다. 짧은 거리에서는 두 스핀이 공간 궤도를 공유하는 해(RKS)와 스핀별 궤도를 독립적으로 푼 해(UKS)의 에너지가 일치했다. RKS 해가 스핀을 반대로 회전시키는 방향에서 불안정해지는 거리는 **1.6178125–1.618125 Å**였다. 그림에서 점선은 그 구간의 중앙이다. 이 거리는 선택한 함수·기저·격자에서의 평균장 분기점이다. [전체 에너지·스핀 표](validation/h2_scan.csv) · [해리 곡선과 스핀 곡률](validation/h2_dissociation.png)

| H–H (Å) | UKS − RKS (mHa) | RKS의 최저 스핀 곡률 (Ha) | 한 H의 절대 Mulliken 스핀 분포 (e) | UKS ⟨S²⟩ |
| ---: | ---: | ---: | ---: | ---: |
| 0.74 | 0.000000 | +0.630895 | 0.0000 | 0.0000 |
| 1.61 | 0.000000 | +0.003066 | 0.0000 | 0.0000 |
| 1.62 | −0.000881 | −0.000843 | 0.0975 | 0.0083 |
| 1.70 | −1.031221 | −0.030229 | 0.5222 | 0.2450 |
| 2.00 | −14.120559 | −0.115384 | 0.8660 | 0.7168 |
| 4.00 | −77.588935 | −0.273903 | 0.9997 | 0.9997 |

모든 분자 계산은 고정 점핵·비상대론·전하 0, α 전자 1개와 β 전자 1개, 적분 격자 level 3이다. UKS의 총 스핀 투영 Mₛ는 0으로 고정된다. 결합을 늘리면 두 원자에 반대 스핀 밀도가 나타나고, 제한을 푼 단일 행렬식에 스핀 오염이 생긴다. 4 Å의 UKS ⟨S²⟩ ≈ 1은 순수 단일항(0)과 순수 삼중항(2) 사이에 있다. 이 스핀 대칭 깨짐은 근사 전자구조 계산의 성질로 해석한다. [Hait·Rettig·Head-Gordon의 H₂ 해리 연구](https://arxiv.org/abs/1810.09083)

자체 RKS는 DIIS로 수렴시켰다. RKS 공간 궤도를 UKS 두 스핀에 복사한 뒤, 정규화한 반대 방향 회전 P = (I, −I)ᵀ/√2에 전체 실수 궤도 Hessian을 투영해 PᵀHP의 최저 고유값을 구했다. 이 고유값이 음수이면 대칭 해에서 에너지를 낮출 수 있다. UKS에는 해당 모드로 초기 궤도를 회전시키고 기존 Newton 최소화기를 적용했다. 시범 계산에서는 1.70 Å의 UKS-DIIS가 수렴 판정을 받고도 대칭 해로 돌아갔다. 낮은 에너지 해는 위 표에 기록했다.

PySCF 2.14.0의 별도 RKS 및 UKS SCF로 각 점을 확인했다. PySCF UKS의 초기 밀도는 각 원자의 1s 원자 궤도에 반대 스핀을 배치해 만들었다. 두 엔진은 같은 AO 적분과 Libxc PBE 및 같은 적분 격자를 사용한다. 분기점 양 끝에서는 PySCF UKS 궤도 기울기를 회전 좌표에서 ±10⁻⁴로 중앙차분해 전체 반대 스핀 곡률 행렬을 다시 구성했다. 양 끝점의 최저 곡률은 자체 계산에서 +7.17×10⁻⁶, −1.145×10⁻⁴ Ha였고 PySCF 중앙차분도 같은 부호였다. [재현 스크립트](../scripts/validate_h2_dissociation.py)

19점에서 두 엔진 사이의 총에너지 최대 차이는 5.64×10⁻¹⁴ Ha, UKS ⟨S²⟩의 최대 차이는 1.23×10⁻⁶였다. 적분 격자를 level 4로 바꿔도 같은 양 끝에서 곡률은 +7.33×10⁻⁶, −1.143×10⁻⁴ Ha였다. 같은 기저·함수로 계산한 고립 H 이중항 두 개의 합은 −0.9972588924 Ha였으며 4 Å의 UKS는 이 값보다 0.00003240 Ha 낮았다. 두 차례 거리 스캔의 CSV가 바이트 단위로 일치했고, WSL에서 `test_stability`와 `test_native_response`의 기존 시험 4건이 통과했다.

```powershell
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_h2_dissociation.py --out outputs/h2_dissociation/recheck
```

스크립트는 새 폴더에 거리별 표, 이분법 이력과 곡률 대조 결과, 그림을 저장한다. 여기서 검증한 스핀 반대 회전 검사는 H₂ 연구 스크립트에 들어 있다. 일반 분자용 `native_stability=check|follow` 옵션은 기존 실수 내부 안정성 검사를 수행한다.

## 질소 기체의 열용량·엔탈피·엔트로피 (1 bar)

QS-20261121-N2-THERMO. 중성 질소 분자의 전자 바닥상태에서 자체 RKS PBE/def2-TZVP로 구조를 최적화하고 조화 진동수를 구했다. 이 구조와 진동수를 고정해 이상기체·강체 회전자·조화 진동자(RRHO)로 298.15 K와 300–2,000 K의 50 K 간격, 총 36개 온도를 계산했다. [전체 표 CSV](validation/n2_thermo.csv) · [온도별 편차 그림](validation/n2_errors.png)

표의 괄호는 계산값에서 NIST Shomate 기준을 뺀 편차다. ΔH는 같은 기체의 H(T)−H(298.15 K)이다. NIST의 100–500 K 계수를 500 K 미만에, 500–2,000 K 계수를 양 끝점 포함 구간에 적용했다. 계수 반올림으로 298.15 K에서 남는 0.0000276 kJ/mol을 모든 NIST 엔탈피에서 빼 두 표의 영점을 맞췄다. [NIST 질소 기체 열물성](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1#Thermo-Gas)

| T (K) | Cp J/(mol·K) | ΔH kJ/mol | S J/(mol·K) |
| ---: | ---: | ---: | ---: |
| 298.15 | 29.113 (−0.011) | 0.000 (0.000) | 191.640 (+0.031) |
| 500 | 29.538 (−0.043) | 5.906 (−0.006) | 206.758 (+0.018) |
| 1,000 | 32.556 (−0.136) | 21.410 (−0.053) | 228.127 (−0.044) |
| 1,500 | 34.629 (−0.218) | 38.263 (−0.142) | 241.766 (−0.114) |
| 2,000 | 35.684 (−0.292) | 55.870 (−0.266) | 251.889 (−0.185) |

계산한 Cp는 36개 온도 모두에서 기준보다 작았다. 2,000 K에서 차이가 0.81%까지 커졌고, 298.15 K부터 누적된 엔탈피 증가량을 0.47% 작게 예측했다. 엔트로피 편차는 650 K에서 양수, 700 K에서 음수였다. 298.15 K에서 별도로 비교한 CODATA 검토값은 191.609 ± 0.004 J/(mol·K)이며, 계산값 191.6399와의 차이 +0.0309는 이 불확도보다 크다. 계산은 질소-14의 질량 14.003074 u, 회전 대칭수 2, 단일항 다중도 1을 사용했다. 비조화 진동, 회전·진동 결합 및 들뜬 전자 상태의 기여를 따로 추정하지 않았다. [NIST의 CODATA 기준과 Shomate 표](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1#Thermo-Gas)

최적화는 geomeTRIC tight, 수치적분 격자 level 4, SCF 에너지 허용오차 10⁻¹⁰ Ha로 수행했다. 평형 결합 길이는 1.102925 Å, 최대 핵 기울기는 2.31×10⁻⁸ Ha/bohr였다. 0.005 bohr 중심차분으로 진동 모드 1개(2353.73 cm⁻¹)를 얻었다. NIST 바닥상태의 분광 상수 rₑ = 1.097685 Å, ωₑ = 2358.57 cm⁻¹와 대조하면 각각 +0.00524 Å, −4.84 cm⁻¹이다. [NIST 이원자 질소 분광 상수](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1000)

선형 분자의 Cp는 7R/2 + Cv,vib로 계산했다. 일정 압력에서 독립적으로 중앙차분한 dH/dT 및 T·dS/dT와의 최대 차이는 1.14×10⁻⁸ J/(mol·K)였다. 300, 1,000, 2,000 K에서 Shomate 구현을 NIST의 반올림 표와 대조했다. 재현 시 두 명령은 각각 새 출력 폴더를 사용한다. [NIST 온도별 표](https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1&Table=on&Type=JANAFG)

변위를 0.0025 bohr로 줄이면 진동수는 2353.6820 cm⁻¹로 0.0444 cm⁻¹ 낮아졌다. 36개 온도에서 두 변위 간 최대 차이는 Cp 0.000106 J/(mol·K), ΔH 0.000129 kJ/mol, S 0.000125 J/(mol·K)였다. 따라서 관측된 2,000 K 편차는 이 변위 변경에 따른 차이보다 크다. WSL에서 기존 열역학 검사 4개(ASE·PySCF 대조 포함)가 통과했다.

```powershell
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_n2_thermo.py --out outputs/n2_thermo/recheck
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_n2_thermo.py --freq-step-bohr 0.0025 --out outputs/n2_thermo/recheck_halfstep
```

## 0.9.8 · 원자 SCF 혼합 방법 비교

QS-20261121-ATOM-MIX. 중성 H–Kr의 구면 평균 바닥 배치 36종에 대해 스핀 비분극·분극 PBE를 각각 계산했다. 두 방법 모두 기본 격자(h = 0.001, r_max = 50 bohr), Slater 초기 밀도, 감쇠 0.3, 동일한 에너지·밀도 수렴 기준을 사용했다. 72쌍 모두 첫 반복 기록이 일치했고, 두 방법 모두 수렴했다.

| 전체 72조건 | 선형 | Anderson |
| --- | ---: | ---: |
| SCF 반복 범위 | 37–58회 | 8–23회 |
| 반복 횟수 합계 | 3,768회 | 1,060회 |
| 반복 횟수 중앙값 | 53.5회 | 15회 |
| SCF 경과 시간 합계 | 296.84초 | 87.35초 |

72조건 모두 반복 수와 실측 시간이 줄었다. 시간 합계의 비는 3.40배다. Anderson 실행의 첫 선형 갱신을 제외하고 음의 밀도 때문에 15회, 불안정한 최소제곱 이력 때문에 4회 선형 혼합으로 복귀했다. 반복 수가 가장 많았던 분극 Cu는 58→23회였다. [전체 72쌍 원시 기록](../outputs/atom_mixing_validation/20261121_h_kr.jsonl)

두 방법 사이의 최대 차이는 총에너지 3.64×10⁻¹² Ha(분극 Kr), 스핀별 방사 밀도의 적분 절대 차이 합 3.11×10⁻⁸ e(분극 C), 점유 오비탈 고윳값 7.08×10⁻⁸ Ha(분극 Zn)였다. 현재 선형 계산의 반복 수는 0.9.7의 WSL 기록과 72조건 모두 같았고 총에너지 차이는 최대 1.37×10⁻¹² Ha였다. 같은 Windows 환경에서 이전 커밋 `9fb2c0a`와 직접 비교한 He 비분극과 B 분극(h = 0.004, LDA 및 PBE)의 반복 이력·에너지 성분·오비탈 값·최종 스핀 밀도 배열은 모두 정확히 일치했다.

Windows 11, Python 3.13.5, NumPy 2.5.4, SciPy 1.18.1에서 수치 라이브러리 스레드를 1개로 고정했다. 두 방법을 예열하고 계산 순서를 바꿔 가며 `run_atom`의 경과 시간만 쟀다. 72조건은 각 방법을 한 번씩 측정했다. 시간은 이 장비와 계산 조건에 대한 값이다. [계산 방법과 안전장치](ATOMIC_PBE_KO.md)

대표 4조건은 별도 실행에서 세 번씩 반복했다. 아래 시간은 각 방법의 중앙값이며, 반복 횟수는 세 번 모두 같았다. 밀도 차이는 두 최종 출력 밀도의 스핀별 적분 절대 차이 합이다. [대표 원자 반복 측정 원본](../outputs/atom_mixing_validation/20261121_pilot.jsonl)

| 원자와 스핀 | 반복 (선형 → Anderson) | 시간 (초, 선형 → Anderson) | 총에너지 절대 차이 (Ha) | 밀도 차이 (e) |
| --- | ---: | ---: | ---: | ---: |
| He 비분극 | 41 → 8 | 0.573 → 0.114 | 1.33×10⁻¹⁵ | 1.92×10⁻⁹ |
| B 분극 | 53 → 12 | 4.580 → 1.035 | 1.07×10⁻¹⁴ | 1.52×10⁻⁸ |
| Fe 분극 | 56 → 16 | 7.819 → 2.285 | 0 | 1.69×10⁻⁸ |
| Kr 비분극 | 54 → 15 | 2.926 → 0.819 | 4.55×10⁻¹³ | 1.91×10⁻⁸ |

```powershell
.\.venv\Scripts\python.exe scripts/validate_atom_mixing.py --out outputs/atom_mixing_validation/recheck_all.jsonl
.\.venv\Scripts\python.exe scripts/validate_atom_mixing.py --representative --repeats 3 --out outputs/atom_mixing_validation/recheck_four.jsonl
```

각 명령은 새 파일에 조건별로 한 줄씩 기록한다. `passed`는 양쪽 수렴, 첫 반복 일치, |ΔE| < 10⁻⁹ Ha, 스핀 밀도 차이 < 2×10⁻⁷ e, Anderson 전자 수 오차 < 10⁻⁹ e, 양쪽 품질 경고 0건을 검사한다.

## 0.9.7 · 구면 평균 원자 PBE

11월 18일 기존 NIST LDA/LSDA 72건을 두 격자에서 다시 계산해 판정 기준을 모두 통과했다. Richardson 외삽값과 NIST 표의 최대 절대 차이는 총에너지 4.99×10⁻⁷, 오비탈 고윳값 5.85×10⁻⁷, 에너지 성분 7.87×10⁻⁷ Ha였다. [전체 재검사 결과](../outputs/atom_validation/20261118_pbe_regression/report.json)

QS-20261117-ATOM-GGA. 자체 PBE로 H–Kr 바닥 배치의 비분극·분극 72건이 기본 로그 격자(h=0.001)에서 37–58회 만에 모두 수렴했다. 최대 전자 수 오차는 1.43×10⁻¹⁴ e, 직접식과 이중 계산 보정식의 최대 차이는 1.82×10⁻¹² Ha였다. H·N·Fe 분극 및 He·Ar·Kr 비분극에서 h=0.002, 0.001, 0.0005의 관측 수렴 차수는 1.999994–2.000002였다. 각 오비탈의 마디 수를 확인했다. [수식·격자·재현 명령](ATOMIC_PBE_KO.md)

표의 Δ는 h=0.002와 0.001에서 구한 방사 격자 Richardson 외삽값에서 독립 기준값을 뺀 값이다. 판정 기준은 |Δ|≤10⁻⁶ Ha이다. 앞의 다섯 원자는 PySCF 2.14.0 RKS와 Libxc 7.0.0, Zn·Kr은 Libxc와 구대칭 Gaussian Galerkin RKS를 사용했다.

| 원자 | 독립 기준 총에너지 (Ha) | Δ (Ha) |
| --- | ---: | ---: |
| He | −2.892934865842 | −1.01×10⁻⁹ |
| Be | −14.629947189590 | −5.26×10⁻⁷ |
| Ne | −128.866427699980 | −4.46×10⁻⁸ |
| Mg | −199.955114950434 | −2.18×10⁻⁷ |
| Ar | −527.346127806298 | −9.68×10⁻⁷ |
| Zn | −1779.182796664558 | −4.66×10⁻⁸ |
| Kr | −2753.416108819394 | −1.16×10⁻⁷ |

열린 껍질 붕소(B, 2p¹)는 독립 PySCF UKS의 `frac_occ_` 실행에서 두 기저 모두 다수 스핀 2p 점유가 (⅓, ⅓, ⅓)이었다. 반지름 0.002–8 bohr의 25개 점·5개 방향에서 매 반복과 최종 스핀 밀도의 최대 방향 차이는 7.21×10⁻¹¹ e/bohr³였다. 두 Gaussian 기저의 에너지 차이는 3.96×10⁻⁸ Ha, Richardson −24.605381482867 Ha와 최종 UKS의 차이는 −1.83×10⁻⁷ Ha였다. `outputs/atomic_pbe_validation/20261118_fractional/fractional/report.json`에 실제 점유·각도별 차이·수렴 판정을 저장했다. 이 비교의 범위는 B의 구면 평균 배치이다.

PySCF RKS는 짝수 템퍼드 지수비 1.4→1.3, 겹침 고유값 절단 10⁻¹⁰, 원자당 600×86 적분점, SCF 에너지 기준 10⁻¹⁰ Ha를 사용했다. 두 기저 모두 수렴했고 기저 변화는 최대 1.16×10⁻⁷ Ha였다. Ar은 허용오차에 가장 가깝다. 독립 구대칭 가우시안 기준으로 다시 계산한 Ar의 |Δ|는 3.1×10⁻⁸ Ha였다.

Zn·Kr은 별도 검증 스크립트에서 s, p, d 가우시안 기저, 약형 GGA 행렬, Gauss–Legendre 부분구간 적분을 사용했다. 1.35→1.25 지수비 변화의 최대치는 두 원소에서 9.13×10⁻⁸ Ha, 70→100 구간 변화는 6.83×10⁻¹² Ha 이하였다. 매우 큰 가상상태 에너지 때문에 첫 겹침 기준 직교화에서 He의 밀도 잔차가 정체했다. 직교화에 운동에너지와 Z² 곱하기 겹침행렬을 더한 양의 계량을 사용해 같은 수렴 조건을 충족했다. 비교 스크립트는 제품의 방사 고유상태, 하트리, XC 구현을 가져오지 않는다.

문헌과의 별도 대조에서는 Ryan 등(2022)의 비상대론 PBE 표 A1에 N −54.5358 Ha, 표 A3에 Zn −1779.1840 Ha 및 Kr −2753.4175 Ha가 실려 있다. 구면 평균 방사 외삽값과의 차이는 각각 약 +0.0000445, +0.00120, +0.00139 Ha다. 이 표는 소수 넷째 자리까지 인쇄되었으며 위의 10⁻⁶ Ha 판정에는 독립 계산값을 사용했다. [원문 표 A1·A3](https://doi.org/10.4236/jamp.2022.1011225)

LDA/LSDA는 시작 커밋 `26381a7` 원본과 같은 Windows 프로세스에서 H–Kr 72건의 에너지·성분·고윳값·SCF 이력·스핀 밀도·퍼텐셜·파동함수가 모두 동일했다. 과거 WSL 보고서와의 총에너지 차이는 최대 3.19×10⁻¹² Ha였다. 관련 WSL unittest 58개 통과, Windows 54개 통과·외부 라이브러리 시험 4개 skip. CLI 설정 재실행과 PBE 미수렴 종료 코드 2도 확인했다.

원시 보고서는 `outputs/atomic_pbe_validation/`의 `20261117_v1/pbe`, `20261117_v2/lda-source`, `20261117_v2/reference`, `20261117_final/compare`, `20261118_gaussian_v2` 아래에 있다. 시험 로그는 `outputs/cell_validation/atom_pbe_regression` 및 `atom_pbe_windows_final`에 보존한다.

11월 18일 최종 WSL 회귀 58개를 다시 통과했고, B·Kr 외부 기준 에너지를 고정하는 시험을 추가한 뒤 영향받는 SCF 시험 4개도 통과했다(고유 시험 59개). Windows에서도 신규 B·Kr 및 수정된 시험을 확인했다(고유 55개 통과, 4개 skip). WSL 최종 로그는 `outputs/cell_validation/atom_pbe_final_20261118`와 `atom_pbe_external_energy_tests`에 있다. Windows CLI의 0.9.7 B PBE 실행은 53회에 수렴했고 에너지 −24.605382659330 Ha, 적분 전자 수 5.000000000000이었다.

## 0.9.6: 원자 LDA/LSDA (로그 방사 격자, H–Kr)와 NIST·PySCF 대조

작업 ID: QS-20261117-ATOM-LSDA, Claude. 시작 HEAD `79856af`. 환경·의존성 변경 없음. Windows Python 3.13.5 / NumPy 2.5.4 / SciPy 1.18.1, WSL Python 3.12.3 / PySCF 2.14.0 / Libxc 7.0.0. 실행 2026-11-17. GUI 시험은 `QUANTUMSTARTER_TK_HIDDEN=1`.

**모델과 가정.** 원자 단위, 고정 점핵, 비상대론, 구면 평균 밀도(부분 껍질은 m 성분 균등 점유 — NIST 표와 같은 관례). LDA는 스핀 비분극, LSDA는 각 껍질을 다수 스핀부터 2l+1개까지 채운다(Hund 최대 스핀). 기본 배치는 NIST 표 머리글(Z = 1–36; Cr `[Ar] 3d5 4s1`, Cu `[Ar] 3d10 4s1`). XC: Slater 교환 + VWN5 상관(Libxc `LDA_X` 1 + `LDA_C_VWN` 7). LSDA 상관은 VWN 스핀 보간 ε_c = ε_P + α_c f(ζ)(1−ζ⁴)/f″(0) + (ε_F − ε_P) f(ζ)ζ⁴(강자성 A = 0.01554535, b = 7.06042, c = 18.0578, x₀ = −0.325; 스핀 강성 A = −1/(6π²), b = 1.13107, c = 13.0045, x₀ = −0.0047584)를 직접 구현했다. 비분극 XC는 0.9.5 `atomic_dft`의 식을 그대로 쓴다.

**수치 방법.** 격자 r = a(eˣ − 1), x_i = i·h, a = 10⁻⁴/Z bohr, 바깥 경계 ≈ 50 bohr(Dirichlet). u = J^½φ(J = dr/dx = r + a)로 바꾸면 −½φ″ + [1/8 + J²(l(l+1)/(2r²) + v)]φ = εJ²φ가 정확히 성립하며(Liouville 변환), 2차 중앙차분 후 C = J⁻¹AJ⁻¹로 대칭 삼중대각화해 SciPy `eigh_tridiagonal`(stebz)로 푼다. 운동에너지에는 변환에서 나온 1/8 항이 들어간다(∫u′²dr/2 = ∫φ′²dx/2 + ∫φ²dx/8). 적분은 x에 대한 사다리꼴(가중 J·h), 하트리는 안쪽 전하 Q(r)/r + 바깥 구껍질 ∫ρ/r′dr′의 누적합이다. SCF는 스핀별 방사 밀도의 선형 혼합(기본 0.3), 잔차는 혼합 전 Σ_σ∫|ρ_out,σ − ρ_in,σ|dr, 총에너지는 고윳값 합의 이중 계산 보정식이며 직접 식(T_s + E_en + E_H + E_xc)과 1e-14 Ha 내에서 같다. 초기 밀도는 Slater 규칙 유효 핵전하의 수소형 오비탈이다.

**발견 — 계단형 행렬의 고유값 정밀도.** C는 원점 근처 성분이 ~1/(h·a)²(예: Z = 18, h = 0.01에서 3e14)라 `eigh_tridiagonal`의 기본 허용오차 eps·‖C‖₁로는 고윳값 오차가 0.07 ~ 10³ Ha였다(scratch 탐색, a = 10⁻⁴/Z ~ 10⁻⁸/Z, Z = 1/18/36, l = 0/1). 이분법(stebz)에 `tol = 4·tiny`를 주면 상대 정밀도로 수렴해 Rayleigh 몫과 상대 ~1e-12로 일치했다. 엔진은 이 설정을 쓰고 ε는 차분 제곱 합 형태의 Rayleigh 몫으로 다시 계산한다(원 고유값과의 차는 `eigensolver_minus_rayleigh_hartree`에 기록, ~1e-9 Ha 이하). 순수 쿨롱(Z = 18, l = 0–2)에서 고윳값 상대 오차 ≈ 0.05·h², h 절반마다 오차 1/4, 비리얼 T = −E, 정규화 1e-12, 마디 수 n − l − 1을 시험한다.

**발견 — 기본 격자 결정.** He LDA에서 h = 0.01/0.005/0.0025/0.00125의 E는 −2.834837176/−2.834836012/−2.834835721/−2.834835648 Ha(차이 비 4.0)이고, 마지막 두 값의 Richardson −2.834835624055는 0.9.5 균일 격자 Richardson과 같으며 PySCF −2.834835624033과 2e-11 차이다. h = 0.0025에서 Ar 오차가 1.5e-4 Ha여서 기본 h를 0.001로 정했다(점 ~11,000–17,000개, 원자당 0.5–10초).

**NIST 기준 데이터.** S. Kotochigova, Z. H. Levine, E. L. Shirley, M. D. Stiles, C. W. Clark, *Atomic Reference Data for Electronic Structure Calculations* (NIST), https://math.nist.gov/DFTdata/atomdata/tables/ptable.html. 2026-11-17에 Z = 1–36 페이지를 내려받아(1초 간격) `scripts/build_nist_atomic_reference.py`가 `<PRE>` 표를 고정 패턴으로 파싱해 `data/nist_atomic_reference.json`(LDA·LSD 열의 Etot·Ekin·Ecoul·Eenuc·Exc와 궤도 고윳값, 페이지 SHA-256)을 만들었다. 숫자는 손으로 옮기지 않았다. notation 페이지에 따라 LSD 고윳값은 다수 스핀이 먼저, 소수 스핀이 다음 줄 LSD 열에 있으며 파서가 열 위치로 확인한다. N LSD의 궤도 고윳값 6개가 NIST와 6자리 일치해 NIST LSD가 같은 VWN 스핀 보간과 껍질별 다수 스핀 우선 점유를 쓴다는 것을 확인했다.

**검증 스크립트.** `wsl --exec .venv-dft/bin/python scripts/validate_atoms.py --label 20261117_h_kr_v3`, exit 0, `outputs/atom_validation/20261117_h_kr_v3/`(report.json, nist_comparison.csv, progress.log). 같은 코드의 첫 두 실행 `20261117_h_kr`, `20261117_h_kr_v2`는 Libxc 대조 항목의 분리 기준만 다르고 NIST·PySCF 값은 비트 단위로 같다(보존). 각 원소·LDA/LSD를 h = 0.002와 0.001로 계산하고 Richardson X(0.001) + [X(0.001) − X(0.002)]/3로 h² 항을 없앴다. 판정 기준(H–Ar): 모두 수렴, 기본 |ΔE| ≤ 1e-4, Richardson |ΔE| ≤ 2e-6, Richardson max |Δε| ≤ 2e-6, PySCF |ΔE| ≤ 1e-6 Ha. 결과: 72건 모두 수렴, 모든 오비탈 마디 수 정상, 판정 모두 통과.

| 구간 | 기본 격자 max \|ΔE_tot\| | Richardson max \|ΔE_tot\| | Richardson max \|Δ성분\| | 기본 max \|Δε\| | Richardson max \|Δε\| |
| --- | --- | --- | --- | --- | --- |
| H–Ar (36건) | 2.4e-05 | 5.0e-07 | 5.1e-07 | 2.8e-06 | 5.5e-07 |
| K–Kr (36건) | 1.6e-04 (상대 5.7e-08) | 4.7e-07 | 8.0e-07 | 1.8e-05 | 5.8e-07 |

NIST 표는 소수 6자리(반올림 폭 ±5e-7 Ha)이므로 Richardson 차이는 반올림 폭과 같은 크기다. 세 번째 격자 h = 0.0005로 본 관측 차수는 N LSD·Ar LDA·Kr LSD 모두 2.000이고, (0.001, 0.0005) Richardson과 NIST의 차이는 +3.4e-7, +8.1e-8, −4.2e-7 Ha다.

원소별(ΔE = 계산 − NIST [Ha]):

| Z | 원소 | 배치 | LDA 기본 ΔE | LDA Richardson ΔE | LSD 기본 ΔE | LSD Richardson ΔE | Richardson 최대 \|Δε\| (LDA / LSD) | 반복 (LDA / LSD) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | H | 1s1 | +4.8e−07 | +4.8e−07 | +2.4e−07 | +2.4e−07 | 4.3e-09 / 2.2e-07 | 36 / 37 |
| 2 | He | 1s2 | +3.6e−07 | +3.8e−07 | +3.6e−07 | +3.8e−07 | 2.8e-07 / 2.8e-07 | 41 / 41 |
| 3 | Li | [He] 2s1 | −3.5e−07 | −1.9e−07 | +1.4e−07 | +3.1e−07 | 3.4e-07 / 4.1e-07 | 47 / 48 |
| 4 | Be | [He] 2s2 | −9.3e−07 | −4.7e−07 | −9.3e−07 | −4.7e−07 | 3.9e-07 / 3.9e-07 | 47 / 47 |
| 5 | B | [He] 2s2 2p1 | −9.4e−07 | −1.0e−07 | −5.1e−07 | +3.3e−07 | 1.7e-07 / 4.1e-07 | 48 / 52 |
| 6 | C | [He] 2s2 2p2 | −8.5e−07 | +4.6e−07 | −9.6e−07 | +3.4e−07 | 2.9e-07 / 4.9e-07 | 47 / 52 |
| 7 | N | [He] 2s2 2p3 | −2.0e−06 | −1.4e−07 | −1.5e−06 | +3.4e−07 | 3.1e-07 / 3.8e-07 | 47 / 51 |
| 8 | O | [He] 2s2 2p4 | −2.3e−06 | +2.0e−07 | −2.3e−06 | +1.9e−07 | 2.2e-07 / 4.9e-07 | 46 / 52 |
| 9 | F | [He] 2s2 2p5 | −3.3e−06 | −1.1e−07 | −2.7e−06 | +4.7e−07 | 3.2e-07 / 4.9e-07 | 46 / 49 |
| 10 | Ne | [He] 2s2 2p6 | −4.2e−06 | −2.7e−07 | −4.2e−06 | −2.7e−07 | 4.3e-07 / 4.3e-07 | 46 / 46 |
| 11 | Na | [Ne] 3s1 | −5.6e−06 | −3.1e−07 | −5.7e−06 | −4.0e−07 | 5.2e-07 / 5.2e-07 | 50 / 52 |
| 12 | Mg | [Ne] 3s2 | −7.4e−06 | −3.2e−07 | −7.4e−06 | −3.2e−07 | 4.3e-07 / 4.3e-07 | 49 / 49 |
| 13 | Al | [Ne] 3s2 3p1 | −9.5e−06 | −4.1e−07 | −9.2e−06 | −8.7e−08 | 4.2e-07 / 5.0e-07 | 51 / 52 |
| 14 | Si | [Ne] 3s2 3p2 | −1.1e−05 | +4.0e−07 | −1.1e−05 | +1.7e−07 | 4.4e-07 / 4.3e-07 | 51 / 52 |
| 15 | P | [Ne] 3s2 3p3 | −1.4e−05 | +3.9e−08 | −1.4e−05 | −2.7e−07 | 2.5e-07 / 5.5e-07 | 51 / 52 |
| 16 | S | [Ne] 3s2 3p4 | −1.7e−05 | +5.0e−07 | −1.7e−05 | +2.4e−07 | 2.8e-07 / 4.6e-07 | 51 / 51 |
| 17 | Cl | [Ne] 3s2 3p5 | −2.1e−05 | −4.9e−07 | −2.1e−05 | −6.1e−08 | 3.6e-07 / 4.0e-07 | 51 / 51 |
| 18 | Ar | [Ne] 3s2 3p6 | −2.4e−05 | +8.1e−08 | −2.4e−05 | +8.1e−08 | 4.8e-07 / 4.8e-07 | 51 / 51 |
| 19 | K | [Ar] 4s1 | −2.9e−05 | +2.9e−07 | −2.9e−05 | −2.2e−07 | 5.2e-07 / 4.5e-07 | 54 / 54 |
| 20 | Ca | [Ar] 4s2 | −3.4e−05 | +3.9e−07 | −3.4e−05 | +3.9e−07 | 5.2e-07 / 5.2e-07 | 54 / 54 |
| 21 | Sc | [Ar] 3d1 4s2 | −3.9e−05 | −3.7e−07 | −3.9e−05 | +4.5e−07 | 4.5e-07 / 5.2e-07 | 54 / 54 |
| 22 | Ti | [Ar] 3d2 4s2 | −4.4e−05 | +1.1e−07 | −4.4e−05 | +4.3e−07 | 4.7e-07 / 4.4e-07 | 54 / 54 |
| 23 | V | [Ar] 3d3 4s2 | −5.0e−05 | −3.2e−07 | −4.9e−05 | +1.4e−07 | 5.2e-07 / 4.6e-07 | 54 / 54 |
| 24 | Cr | [Ar] 3d5 4s1 | −5.4e−05 | +3.6e−08 | −5.4e−05 | −4.1e−08 | 4.8e-07 / 5.0e-07 | 53 / 54 |
| 25 | Mn | [Ar] 3d5 4s2 | −6.1e−05 | −4.5e−07 | −6.0e−05 | +1.3e−07 | 5.3e-07 / 4.0e-07 | 54 / 55 |
| 26 | Fe | [Ar] 3d6 4s2 | −6.7e−05 | +1.6e−07 | −6.6e−05 | +4.2e−07 | 4.6e-07 / 5.5e-07 | 54 / 56 |
| 27 | Co | [Ar] 3d7 4s2 | −7.3e−05 | −2.1e−07 | −7.3e−05 | −1.5e−07 | 4.9e-07 / 4.9e-07 | 53 / 56 |
| 28 | Ni | [Ar] 3d8 4s2 | −7.9e−05 | +3.0e−07 | −7.9e−05 | +3.0e−07 | 5.3e-07 / 5.0e-07 | 53 / 56 |
| 29 | Cu | [Ar] 3d10 4s1 | −8.5e−05 | +1.3e−07 | −8.5e−05 | +2.6e−07 | 3.8e-07 / 4.1e-07 | 52 / 57 |
| 30 | Zn | [Ar] 3d10 4s2 | −9.3e−05 | +3.2e−07 | −9.3e−05 | +3.2e−07 | 4.2e-07 / 4.2e-07 | 53 / 53 |
| 31 | Ga | [Ar] 3d10 4s2 4p1 | −1.0e−04 | +1.2e−07 | −1.0e−04 | −1.8e−07 | 4.3e-07 / 5.0e-07 | 53 / 53 |
| 32 | Ge | [Ar] 3d10 4s2 4p2 | −1.1e−04 | +5.1e−08 | −1.1e−04 | +2.8e−07 | 4.4e-07 / 5.0e-07 | 52 / 53 |
| 33 | As | [Ar] 3d10 4s2 4p3 | −1.2e−04 | −4.7e−07 | −1.2e−04 | +2.1e−07 | 5.3e-07 / 5.0e-07 | 52 / 53 |
| 34 | Se | [Ar] 3d10 4s2 4p4 | −1.3e−04 | −1.8e−07 | −1.3e−04 | +2.4e−07 | 5.3e-07 / 5.5e-07 | 52 / 53 |
| 35 | Br | [Ar] 3d10 4s2 4p5 | −1.4e−04 | +1.3e−07 | −1.4e−04 | −6.8e−08 | 5.8e-07 / 5.0e-07 | 52 / 52 |
| 36 | Kr | [Ar] 3d10 4s2 4p6 | −1.6e−04 | −4.2e−07 | −1.6e−04 | −4.2e−07 | 4.5e-07 / 4.5e-07 | 52 / 52 |

**독립 기준: PySCF RKS (닫힌 껍질).** `dft.RKS`, xc `LDA_X,LDA_C_VWN`, 비축약 짝수 템퍼드 s(0.02–10⁹) + p(0.02–10⁶, He·Be 제외) 기저, β = 1.8과 1.6, Treutler 300 × Lebedev 50, 가지치기 없음, `small_rho_cutoff = 0`, `conv_tol = 1e-10`(1e-12에서는 거의 선형 종속으로 수렴 판정이 자주 실패해 탐색 후 완화). 기준은 β = 1.6.

| 원자 | PySCF E [Ha] | β 1.8→1.6 변화 | 자체 Richardson − PySCF | 자체 기본 − PySCF | NIST |
| --- | --- | --- | --- | --- | --- |
| He | −2.834835624055 | −2.1e-11 | −2.8e-13 | −1.6e-08 | −2.834836 |
| Be | −14.447209471203 | −1.1e-10 | −2.8e-09 | −4.6e-07 | −14.447209 |
| Ne | −128.233481269076 | −4.3e-09 | −1.4e-10 | −4.0e-06 | −128.233481 |
| Mg | −199.139406294586 | −5.6e-08 | −2.1e-08 | −7.1e-06 | −199.139406 |
| Ar | −525.946194875983 | −4.3e-07 | −4.3e-08 | −2.4e-05 | −525.946195 |

Ar의 자체–PySCF 차이(4.3e-8)는 PySCF 쪽 기저 변화(4.3e-7)보다 작아, 차이의 크기는 PySCF 기저 완전성으로 제한된다고 본다.

**같은 함수 확인(Libxc 7.0.0).** 시험(`test_radial_atom`)은 무작위 스핀 밀도 200점(총밀도 1e-6–1e2, ζ ∈ [−1, 1])에서 교환·상관 각각과 합의 ε, v↑, v↓를 상대 1e-11로 비교하고, ζ = 0에서 비분극 식과 1e-13, v_σ = ∂(nε)/∂n_σ 중앙차분과 1e-6을 확인한다. 수렴한 LSD 밀도(검증 v3)에서는 N·Fe의 최대 상대차가 ε_xc 1.3e-14, v 1.4e-14(N 소수 스핀 1.0e-13)다. Cr(3d5 4s1)은 소수 스핀에 3d·4s 전자가 없어 바깥 꼬리에서 소수 스핀 밀도가 1e-27 bohr⁻³까지 작아진다. 두 스핀 밀도가 모두 1e-15 이상인 점에서는 ε 1.1e-15, v 1.4e-11로 같지만, 소수 스핀 밀도가 1e-15 미만인 561점(전자 1.1e-3개)에서는 소수 스핀 퍼텐셜이 최대 1.6% 다르다. 조사 결과 n↑ = 1e-6에서 Libxc의 n↓ = 1e-16 값이 자체 식의 n↓ = 1e-15 값과 10자리까지 같았다. 즉 Libxc 7은 1e-15 미만 스핀 밀도를 1e-15로 바꿔 계산한다. 자체 값은 한쪽 차분 극한과 일치하며(최악 점 자체 −7.0919e-4, 차분 −7.0920e-4, Libxc −6.980e-4), Libxc의 ε를 대신 써도 E_xc는 1.8e-15 Ha만 바뀐다. 처음에는 ζ 문턱값으로 추정해 1−|ζ| ≤ 1e-12를 분리했으나(검증 v2) 그 기준으로는 2.9e-3 차이가 남아 밀도 문턱값으로 정정했다(검증 v3, 시험도 같은 기준).

**발견 — 버그 수정.** 첫 CLI 실행에서 에너지가 `numpy.float64`로 남아 수렴 판정이 `numpy.bool_`이 되어 results.json 저장이 TypeError로 실패했다. 에너지와 판정값을 Python float/bool로 바꿨다.

**독립성의 범위.** NIST 값은 별도 원자 코드의 결과(같은 VWN 정의)이고, PySCF는 가우스 기저·해석적 적분·별도 XC 격자를 쓴다. 우리 엔진은 런타임에 PySCF·Libxc를 부르지 않는다. 세 결과가 공유하는 것은 XC 함수의 정의뿐이다.

**한계.** 비상대론·점핵·구면 평균만(다중항 구분, 비구면 밀도, 상대론, GGA·혼성, 음이온 없음). 선형 혼합만 있어 반복이 36–57회다. 기본 격자 오차는 |E|에 비례(상대 ~5e-8)하므로 무거운 원자에서 1e-6 Ha 수준이 필요하면 h를 줄이거나 두 격자로 외삽한다. Libxc 밀도 문턱값 아래의 꼬리에서 자체 퍼텐셜과 Libxc가 다른 것은 에너지 영향이 없지만, 두 구현의 정의가 그 영역에서 다르다는 점은 기록해 둔다.

**자동 검사.** 신규 `test_radial_atom` 18개: 배치 해석·오류, 기본 배치 = NIST 머리글(36개), Hund 점유, LSDA 식(ζ = 0 극한·미분 일관성·Libxc), 수소형 고윳값 2차 수렴·비리얼·정규화·마디 수, 하트리 해석해, He = PySCF 기준(1e-7)·He LSDA = LDA, N LSD = NIST(총에너지·네 성분 1e-5, 고윳값 2e-6), N Richardson ≤ 2e-6, Ne·Ar LDA = NIST, 전자 수·스핀 모멘트·이중 계산 식, 반복 상한 3 → not_converged 저장, 저장 파일·설정 재로드, CLI 수렴(exit 0)·미수렴(exit 2)·settings.json 재실행 동일 에너지·옵션 충돌 거부. Windows 전체(`QUANTUMSTARTER_TK_HIDDEN=1`, 45.8초): 192개 OK(skip 85 = PySCF·Libxc·실계산 GUI opt-in 필요분). WSL 전체(686초): 172개 중 169 통과, 오류 3개는 `.venv-dft`에 tkinter가 없어 GUI 모듈 3개를 불러오지 못한 기존 환경 제약(0.9.5와 동일). 이 두 실행 뒤 Libxc 시험의 분리 기준만 고쳤으므로 `test_radial_atom test_atomic_dft test_solver test_workflow` 48개를 다시 실행했다: Windows OK(skip 2), WSL OK. 로그 `outputs/atom_validation/windows_all_hidden.log`, `wsl_all.log`.

## 0.9.5: 중성 He 방사 격자 Kohn–Sham LDA-SCF

작업 ID: QS-20261112-HE-LDA, Claude. 시작 HEAD `b0bf6c5`. 환경·의존성 변경 없음(새 패키지 없음). Windows Python 3.13.5 / NumPy 2.5.4 / SciPy 1.18.1, WSL Python 3.12.3 / PySCF 2.14.0 / Libxc 7.0.0 / NumPy 2.5.4 / SciPy 1.18.1. 실행 2026-11-12. GUI 시험은 `QUANTUMSTARTER_TK_HIDDEN=1`.

**모델과 가정.** 원자 단위(Hartree, bohr). 고정 점핵 Z = 2, 비상대론, 전자 2개 1s 이중 점유, 스핀 비분극(제한 닫힌 껍질), 구대칭. u(r) = rR(r), n(r) = 2|u|²/(4πr²), N = ∫n 4πr² dr. 교환·상관: Slater 교환(Libxc `LDA_X`, ID 1, α = 2/3) + VWN5 상관(Libxc `LDA_C_VWN`, ID 7, Ceperley–Alder 상자성 맞춤; A = 0.0310907, b = 3.72744, c = 12.9352, x₀ = −0.10498). PySCF `parse_xc('LDA_X,LDA_C_VWN')` = (1, 7) 확인. 결과 `xc`에 이름·ID·참고문헌을 저장.

**수치 방법.** `solver.py`와 같은 균일 격자 r_i = i·dr, dr = r_max/(N+1), Dirichlet u(0) = u(r_max) = 0, 2차 중앙차분 삼중대각 해밀토니안(SciPy `eigh_tridiagonal`; 0.9.5에서 `radial_grid`·`lowest_radial_states`로 분리, 수소 결과 4개 설정 비트 단위 동일). 하트리 v_H(r) = Q(r)/r + ∫ᵣ ρ(r′)/r′ dr′(ρ = 4πr²n), ρ(0) = ρ(r_max) = 0인 사다리꼴 누적합 — 이산 핵이 dr/max(r_i, r_j)로 대칭이어서 v_H가 이산 E_H의 정확한 미분이다. LDA 식은 닫힌 형태로 직접 계산(n ≤ 1e-30 bohr⁻³에서 0). SCF: 수소형 1s(ζ = 27/16)를 격자에서 정규화한 초기 밀도 → v_eff = −2/r + v_H[n_in] + v_xc[n_in] → 최저 s 상태 → n_out → 혼합 **전** 잔차 ∫|n_out − n_in| 4πr² dr → 선형 혼합. 총에너지는 출력 오비탈의 Kohn–Sham 함수를 고윳값 합과 이중 계산 보정으로 쓴 E = 2ε − ∫(v_H[n_in] + v_xc[n_in]) n_out + E_H[n_out] + E_xc[n_out](수렴 시 2ε − E_H − ∫v_xc n + E_xc). |ΔE|와 잔차가 모두 기준 미만이어야 수렴이다.

**발견 — 고유값 해법의 절대 정밀도.** 첫 구현은 SciPy 고유값을 그대로 ε로 써서 수렴 직전 ΔE가 ±2e-10 Ha에서 요동했다. `eigh_tridiagonal`(stebz) 기본 허용오차가 eps·‖H‖₁ ≈ eps·2/dr²(dr = 0.001에서 4.4e-10 Ha)이기 때문이다. 조치: ε를 고유벡터의 Rayleigh 몫으로 계산하고 운동에너지는 3점 차분과 같은 값인 (f/2)Σ(u_{i+1} − u_i)²/dr로 계산(자릿수 상쇄 없음). 이후 ΔE 잡음은 ~1e-15 Ha, 이중 계산 식과 직접 식(T_s + E_en + E_H + E_xc)의 차이는 7.5e-15 Ha. 원래 고유값과의 차이(~1e-10)는 반복 이력 `eigensolver_minus_rayleigh_hartree`에 남긴다. 같은 이유로 밀도 잔차는 약 0.4·eps·2/dr²(dr = 0.001에서 1.7e-10 e) 아래로 내려가지 않는다(같은 행렬을 다시 풀면 결과 동일 — SCF 문제가 아님). 기본 기준 1e-7 e는 이보다 약 600배 크고, 에너지 오차는 잔차의 제곱에 비례한다(잔차 1e-3에서 5e-11 Ha). 기준이 5·eps·2/dr² 미만이면 경고한다.

**기본 실행(Windows CLI).** `.\.venv\Scripts\python.exe main.py --model he-lda --no-show`, exit 0, version 0.9.5: `outputs/he_lda/20261112_204858_281986_bqah964y/`. 22회 수렴(|ΔE| 1.3e-15 Ha, 잔차 6.2e-8 e), **E = −2.834832194587 Ha(−77.13971379 eV), ε₁ₛ = −0.570423731059 Ha(−15.52202047 eV), N = 2.000000000000**, 0.35초. 성분 [Ha]: T_s 2.767912, E_en −6.625549, E_H 1.996118, E_x −0.861846, E_c −0.111467. 같은 설정의 WSL 결과와 E 4e-15, ε 2e-13 Ha 차(BLAS/LAPACK 빌드 차). `--max-iterations 3`: exit 2, `not_converged`, 마지막 ΔE −5.97e-5 Ha·잔차 3.89e-2 e 저장. 혼합 비율별 반복 수(기준 1e-11 Ha/1e-9 e): 0.3 → 54, 0.5 → 28, 0.7 → 17, 0.8 → 13, 0.9 → 17, 1.0 → 24.

**검증 스크립트.** `wsl --exec .venv-dft/bin/python scripts/validate_he_lda.py --label 20261112_he_lda --compare-result outputs/he_lda/20261112_204858_281986_bqah964y/results.json`, exit 0, `outputs/he_lda_validation/20261112_he_lda/report.json`(기준값 하드코딩 없음). 수열은 |ΔE| 1e-12 Ha·잔차 1e-7 e 기준으로 모두 수렴. 첫 실행 `20261112_final`은 PySCF의 deprecated `remove_linear_dep_`를 썼고 숫자는 같다(보존).

격자 간격 축소(r_max = 30 bohr 고정, PySCF 기준값 대비):

| dr [bohr] | N | 반복 | E [Ha] | E − PySCF [Ha] | 관측 차수 | Richardson E [Ha] | ε₁ₛ [Ha] |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.004 | 7499 | 22 | −2.834780754388 | +5.487e-05 | | | −0.5704089377 |
| 0.002 | 14999 | 22 | −2.834821906274 | +1.372e-05 | | −2.834835623569 | −0.5704207723 |
| 0.001 | 29999 | 22 | −2.834832194587 | +3.429e-06 | 2.0000 | −2.834835624025 | −0.5704237311 |
| 0.0005 | 59999 | 22 | −2.834834766687 | +8.573e-07 | 2.0000 | −2.834835624053 | −0.5704244707 |
| 0.00025 | 119999 | 22 | −2.834835409713 | +2.143e-07 | 2.0000 | −2.834835624055 | −0.5704246559 |
| 0.000125 | 239999 | 22 | −2.834835570470 | +5.356e-08 | 2.0000 | −2.834835624055 | −0.5704247030 |

이산화 오차는 ≈ 3.43·dr² Ha로 2차이며, Richardson 외삽 E(h/2) + [E(h/2) − E(h)]/3은 1e-12 Ha 내에서 일정하다.

계산영역 확대(dr = 0.001 bohr 고정):

| r_max [bohr] | N | E [Ha] | E − E(r_max = 40) [Ha] | ε₁ₛ [Ha] |
| --- | --- | --- | --- | --- |
| 5 | 4999 | −2.8342440304130 | +5.88e-04 | −0.568339930128 |
| 7.5 | 7499 | −2.8348292446606 | +2.95e-06 | −0.570405258061 |
| 10 | 9999 | −2.8348321803405 | +1.42e-08 | −0.570423600700 |
| 12.5 | 12499 | −2.8348321945187 | +6.84e-11 | −0.570423730257 |
| 15 | 14999 | −2.8348321945868 | +3.30e-13 | −0.570423731080 |
| 20 / 30 / 40 | 19999 / 29999 / 39999 | −2.8348321945871 | ≤ 1e-15 | −0.57042373106 ~ 09 |

LDA 퍼텐셜은 지수적으로 0이 되므로 밀도 꼬리는 n ~ e^(−2.14 r)이고, r_max ≥ 15 bohr에서 계산영역 오차는 1e-12 Ha 미만이다. 기본 r_max = 30은 여유값이다.

**같은 함수 확인.** 수렴한 He 밀도 9043점(1e-10 ~ 3.51 bohr⁻³)에서 자체 ε_xc·v_xc와 Libxc 7.0.0 `eval_xc('LDA_X,LDA_C_VWN')`의 최대 상대차는 1.3e-14 / 1.6e-14. 시험(`test_atomic_dft`)에서는 1e-10 ~ 1e4 bohr⁻³ 141점에서 교환·상관 각각과 합을 상대 1e-12로 비교한다. 하트리는 한 전자 수소형 밀도의 해석해 v_H = 1/r − (ζ + 1/r)e^(−2ζr), E_H = 5ζ/16와 비교해 dr = 0.001에서 최대 1.6e-6, E_H 1.0e-7 오차이고 dr 절반마다 1/4로 줄며, 바깥에서 r·v_H = 1(Gauss 법칙)을 1e-12로 만족한다.

**독립 기준: PySCF RKS.** `dft.RKS`, xc `LDA_X,LDA_C_VWN`, 짝수 템퍼드 비축약 s 기저(α_min, β, α_max [bohr⁻²]), Treutler 방사 300점 × Lebedev 50점, 가지치기 없음, `small_rho_cutoff = 0`, `conv_tol = 1e-12`. s 기저만으로 충분한 이유: 구대칭 닫힌 껍질에서 1s 오비탈은 s 대칭이고, p 함수 8개를 더해도 변화가 1e-12 Ha 이하였다.

| α_min, β, α_max | 함수 수 | 최소 겹침 고윳값 | SCF 수렴 | E [Ha] | E − 기준 [Ha] | HOMO [Ha] |
| --- | --- | --- | --- | --- | --- | --- |
| 0.05, 2.5, 1e8 | 25 | 1.9e-03 | 예 | −2.834834054111 | +1.57e-06 | −0.5704255873 |
| 0.05, 2.0, 1e8 | 32 | 8.9e-05 | 예 | −2.834835621394 | +2.64e-09 | −0.5704247290 |
| 0.02, 2.0, 1e8 | 34 | 8.8e-05 | 예 | −2.834835621440 | +2.59e-09 | −0.5704247244 |
| **0.02, 1.8, 1e8 (기준)** | 39 | 8.9e-06 | 예 | **−2.834835624033** | 0 | −0.5704247247 |
| 0.02, 1.6, 1e9 | 54 | 1.8e-07 | 아니오 | −2.834835624055 | −2.1e-11 | −0.5704247679 |

기준 기저에서 방사 75/150/600/1000점, 각도 302점, p 함수 추가에 따른 변화는 최대 6.8e-14 Ha이고, 직전 수렴 기저와의 차이는 2.6e-9 Ha다. β = 1.6은 거의 선형 종속이라 PySCF가 수렴 판정을 내리지 못해 기준에서 뺐다(에너지는 Richardson 값과 같은 −2.834835624055). 참고로 표준 축약 기저는 cc-pVTZ +7.6e-4, cc-pV5Z +1.08e-4, aug-cc-pV5Z +5.0e-5 Ha로, 핵 첨점 때문에 1e-4 Ha 목표에는 큰 지수까지 포함한 기저가 필요했다.

| QuantumStarter | E [Ha] | E − PySCF [Ha] | ε₁ₛ − PySCF HOMO [Ha] |
| --- | --- | --- | --- |
| 기본 (dr = 0.001, r_max = 30) | −2.834832194587 | **+3.43e-06** (목표 ≤ 1e-4 충족) | +9.9e-07 |
| 가장 촘촘한 격자 (dr = 0.000125) | −2.834835570470 | +5.36e-08 | +2.2e-08 |
| Richardson (마지막 두 격자) | −2.834835624055 | −2.2e-11 | +6.0e-09 |
| Windows CLI 저장 결과 | −2.834832194587 | +3.43e-06 | +9.9e-07 |

**문헌값 대조(참고).** NIST Atomic Reference Data for Electronic Structure Calculations(Kotochigova 외)의 He LDA(Vosko–Wilk–Nusair) 값은 Etot = −2.834836, ε₁ₛ = −0.570425, Ekin 2.767922, Eenuc −6.625564, Ecoul 1.996120, Exc −0.973314 Ha다(https://math.nist.gov/DFTdata/atomdata/tables/02He.html, 2026-11-12 열람; 단위는 notation 페이지의 "All data in the tables is given in Hartree energy units"). Richardson 외삽값과 PySCF 기준값은 소수 6자리에서 이와 같다. 기본 격자의 성분은 dr² 오차만큼(≤ 1.5e-5 Ha) 다르다.

**독립성의 범위.** PySCF와 공유하는 것은 XC 함수의 정의(같은 Libxc 식)뿐이다. 우리 엔진은 런타임에 Libxc를 호출하지 않으며, 기저(실공간 격자 vs 가우스), 운동·쿨롱 적분(유한차분·방사 적분 vs 해석적 가우스 적분), XC 적분 격자, 고유값 해법, SCF가 모두 다르다. 함수 정의 자체가 틀렸을 가능성은 NIST 값(별도 원자 코드)과의 일치로 보완된다.

**한계.** 중성 He 1s²만 지원한다(다른 원소·이온, 열린 껍질, 스핀 분극, GGA·혼성 없음). 균일 격자 2차 차분이라 기본 격자의 이산화 오차가 3.4e-6 Ha이고, 격자를 아주 촘촘히 하면 고유값 해법의 부동소수점 하한(eps·2/dr²)이 커지므로 grid_points 상한을 1,000,000으로 두었다. 혼합은 선형뿐이다. LDA의 ε₁ₛ(−0.570 Ha)는 자기 상호작용 오차 때문에 실험 이온화 에너지(24.587 eV ≈ 0.904 Ha)보다 크게 얕다 — 모델의 한계이며 구현 오류가 아니다. 상대론·핵 질량·유한 핵 효과는 포함하지 않는다.

**자동 검사.** 신규 `test_atomic_dft` 18개: LDA 닫힌 식·v = d(nε)/dn·n = 0 처리·Libxc 대조, 하트리 해석해·2차 수렴·안팎 전하, 기본 SCF 수렴(두 기준·변분적 단조 감소), 전자 수 2(1e-10), 이중 계산 식 = 직접 식, PySCF 기준 대비 1e-4 Ha 이내, 격자 2차 수렴과 Richardson ≤ 1e-8 Ha, 계산영역 수렴, 반복 상한 3 → not_converged·마지막 잔차 저장, 저장 파일·설정 재로드·잘못된 입력, CLI 수렴(exit 0)·미수렴(exit 2)·settings.json 재실행 동일 에너지·분자 옵션 거부, `--hydrogen-check` CLI 결과 불변. Windows 전체(`QUANTUMSTARTER_TK_HIDDEN=1`, 22.7초): 174개 OK, skip 84(PySCF·실계산 GUI opt-in 필요분); 수소 GUI `test_gui`는 이번에 추가한 숨김 모드로 2개 통과(화면 픽셀 비교 2줄만 숨김 모드에서 생략). WSL 전체(DISPLAY 해제, 523초): 154개 중 151 통과, 오류 3개는 `.venv-dft`에 tkinter가 없어 GUI 시험 모듈 3개(test_gui·test_molecular_gui·test_periodic_gui)를 불러오지 못한 것 — 바뀌지 않은 import 줄에서 나는 기존 환경 제약이며 이 모듈들은 Windows에서 통과했다. `main.main()`을 쓰는 test_periodic_cell 10개 통과. 로그 `outputs/he_lda_validation/windows_all_hidden.log`, `wsl_all.log`.

## 0.9.4: 분자 IRC와 열역학 독립 교차검증

작업 ID: QS-20261029-NATIVE-DFT 단계 13, Claude. 시작 HEAD `2dd6c8f`. 환경/의존성 설치 변경 없음. WSL Python 3.12.3 / PySCF 2.14.0 / geomeTRIC 1.1.1 / ASE 3.29.0 / NumPy 2.5.4, Windows Python 3.13.5 / ASE 3.29.0. 실행 2026-11-01. GUI 시험은 `QUANTUMSTARTER_TK_HIDDEN=1`(창 표시 없음).

**수치 방법.** `calculation=irc`: geomeTRIC 1.1.1 내장 IRC(Gonzalez–Schlegel 1990 질량가중 2단계, `irc_direction=both`, 끝 근처 최소화 전환). 시작 Hessian은 우리 해석적 힘의 유한차분(6N회). 경로는 geomeTRIC이 받아들인 스텝만(거부 스텝 제외, optimize.py:981 확인). 경로 길이 s는 연속 구조를 질량가중 Kabsch로 겹친 뒤의 거리 합 [bohr·amu^½], 질량은 geomeTRIC이 경로를 정의할 때 쓴 PeriodicTable(표준 원자량). 끝점 판정은 원자별 기울기 노름의 RMS·최대(geomeTRIC calcGradNorm과 같은 정의)를 `opt_convergence` 기준과 비교한다.

**발견 1 — geomeTRIC customengine IRC 결함.** run_optimizer는 IRC일 때 단일 조각 판정에 `M.molecules`를 좌표계 생성 전에 읽는다. 파일로 읽은 분자는 토폴로지가 있지만 customengine 분자는 없어 AttributeError로 중단됐다(재현: 수정 전 첫 NH₃ 실행). 조치: 분자에 `build_topology()`를 먼저 호출(파일 입력 경로와 같은 거리 기반 결합).

**발견 2 — 한 방향의 단계 한도를 geomeTRIC이 성공으로 보고.** HCN⇌HNC 첫 실행(`outputs/irc_validation/hcn_hnc_pbe_svp_20261101/`): HNC 쪽(선형 끝점) 방향이 정상 IRC 23스텝 후 기울기 ~1e-3 Ha/bohr에서 스텝당 ~1e-6 Å로 77스텝 정체, 방향당 한도 100에 도달하자 geomeTRIC은 조용히 반대 방향으로 넘어가 마지막에 "Converged"를 기록했다. 우리 끝점별 검사가 HNC 끝(최대 기울기 9.9e-4 > tight 4.5e-5)을 미수렴으로 판정했다. 조치: 기준 미달 끝점은 그 IRC 구조에서 일반 최소화(`endpoint_<방향>/`)로 마무리하고 IRC 끝점(`irc_end`)·최소화 결과(`relaxation`: 단계, 이동 RMSD, ΔE)를 함께 기록. 이동 RMSD > 0.1 Å이면 경고. 마지막 방향에서 한도에 걸리면 geomeTRIC이 예외로 경로 객체를 잃으므로, 매 평가마다 쓰는 `geometric_irc.xyz`(Å 10자리)를 기록된 계산(전정밀도 에너지·기울기)과 좌표 일치(1e-7 Å)로 맞춰 경로를 복구한다. 전체 수렴은 geomeTRIC 판정이 아니라 두 끝점의 기준 충족으로 정한다(`geometric_converged`는 별도 기록).

**NH₃ 반전 (대칭 재배열, CLI).** 자체 PBE/STO-3G, grid 3, default 기준, TS는 같은 수준 transition_state 결과(`examples/nh3_inversion_ts.xyz`). Windows CLI→WSL `main.py --xyz examples/nh3_inversion_ts.xyz --basis sto-3g --xc pbe --engine native --calculation irc`, exit 0, version 0.9.4: `outputs/irc_validation/cli/20261101_144237_018561_debc906r/`.

| 항목 | 값 |
| --- | --- |
| TS 시작 Hessian 허수 진동 | 1개, −880.60 cm⁻¹ |
| 기울기 계산 / Hessian / 방향별 | 63 / 24 / forward 20·backward 19 |
| 경로 점 | 33 (s = ±1.625 bohr·amu^½) |
| 끝점 에너지 [Ha] | −55.7169351915 / −55.7169356229 (차 4.3e-7 Ha = 0.0011 kJ/mol) |
| 끝점 장벽 [kJ/mol] | 54.985 / 54.986 |
| 끝점 최대 원자 기울기 [Ha/bohr] | 1.42e-4 / 1.36e-4 (기준 4.5e-4) |
| TS에서 오르막 최대 [Ha] | 0 / 0 |
| 두 끝점 RMSD 고유 회전 / 거울상 허용 [Å] | 0.443 / 7.6e-5 (서로 거울상) |

test_irc는 같은 TS에서 독립 최적화한 NH₃와 끝점 에너지(2e-6 Ha)·구조(거울상 허용 RMSD 5e-3 Å) 일치, 방향당 12단계 한도에서 경로 파일 복구 + 양 끝 최소화 후 일치, 최소점 입력의 명시적 거부, geomeTRIC 부재 오류를 검사한다.

**HCN ⇌ HNC (서로 다른 끝점, 연결 확인).** `scripts/validate_irc.py --label hcn_hnc_pbe_svp_20261101_relax`: 자체 PBE/def2-SVP, grid 3, tight 기준(TS 탐색은 default). HCN·HNC 독립 최적화 → 굽은 초기 구조에서 P-RFO TS → TS 진동수 → IRC → 각 끝점을 독립 최적화 구조와 RMSD·에너지로 대조 → 끝점 진동수. 결과 `outputs/irc_validation/hcn_hnc_pbe_svp_20261101_relax/report.json`.

| 항목 | 값 |
| --- | --- |
| TS 허수 진동 (진동수 계산) | 1개, −1077.13 cm⁻¹ |
| IRC | 수렴(끝점 기준), 181.6초, 방향별 계산 forward 103(한도)·backward 34 |
| forward 끝 → HNC | 최소화 4단계(IRC 끝에서 이동 RMSD 4.2e-4 Å, ΔE −6.4e-7 Ha). 독립 HNC와 RMSD 1.29e-4 Å, ΔE 2.0e-8 Ha. 최소화 전 IRC 끝도 HNC로 판별(RMSD 3.9e-4 Å) |
| backward 끝 → HCN | 최소화 불필요. 독립 HCN과 RMSD 5.9e-5 Å, ΔE −5.4e-8 Ha |
| 끝점 진동수 | 둘 다 허수 없음(minimum) |
| 장벽 HCN→TS / HNC→TS | 194.08 / 132.39 kJ/mol |
| 이성질화 HCN→HNC | +61.68 kJ/mol (0.6393 eV) |

에너지 값은 이 계산 수준(PBE/def2-SVP)의 결과이며 문헌·실험값과 비교하지 않았다(기준값 출처 미확인). 검증 대상은 "TS의 IRC가 독립적으로 최적화한 두 이성질체를 잇는다"는 연결과 구현 일관성이다.

**열역학 독립 교차검증.** 같은 진동수·질량·구조·대칭수·다중도를 ASE 3.29 IdealGasThermo에 넣어 우리 `ideal_gas_thermo`와 비교: H₂O(비선형, σ2, 단일항), O₂(선형, σ2, 삼중항), CO₂(선형, σ2), OH(선형, σ1, 이중항), O 원자(삼중항) × (298.15 K/101325 Pa, 1500 K/1e5 Pa, 3000 K/2e6 Pa). 최대 상대차 H 3.39e-7, S 4.09e-7, G는 (|H|+T·S) 대비 3.05e-7. ASE kB(CODATA 2014)와 우리(scipy CODATA) kB의 상대차 −3.39e-7과 같은 크기이므로 상수 판 차이로 설명된다. 첫 시험은 G에 상대 기준 2e-6을 그대로 적용해 물 298 K에서 1.93e-7 eV 차로 실패했다. G=H−TS가 큰 두 항의 차이여서 기준 크기를 |H|+T·S로 바꿨고, 허용 상대값 2e-6은 그대로 두었다. 참고: 물 S(298.15 K, 101325 Pa) = 188.56 J/mol/K.

**자동 검사.** WSL 최종 코드: test_irc/frequencies/optimize/constraints/neb/molecular/energy_analysis/basis_overrides/gradient/dispersion/native/stability 83개 통과(166.8초, `wsl_final.log`). 중간 확인: test_irc 6개(58.3초), 리팩터 직후 회귀 56개(74.4초). Windows: 숨김 실계산 GUI IRC·최적화·진동수 3개(72.8초), 빠른 시험 46개 중 19 통과·27 skip(PySCF/실계산 필요분은 WSL에서 확인). 로그 `outputs/irc_validation/*.log`.

**한계.** 분자 전용(주기계 IRC 없음), 제약 불가, 단일 TS에서 양방향만(재시작 없음). 선형 끝점 근처의 geomeTRIC IRC 정체는 끝점 최소화로 보완하지만 그만큼 계산이 낭비된다(HCN⇌HNC에서 77회). 끝점 최소화는 최급강하 경로의 연장이 아니므로 이동 RMSD를 함께 본다. IRC는 반응 연결의 필요조건 확인이며, 끝점이 의도한 반응물·생성물인지는 사용자가 비교해야 한다. 경로 질량(표준 원자량)과 진동 질량(최빈 동위원소)이 다르다.

## 0.9.3: bulk 셀 완화·FFT mesh 수렴·주기계 GUI

작업 ID: QS-20261029-NATIVE-DFT 단계 12, ChatGPT Codex. 시작 HEAD `2bec7b3`. 환경/의존성 설치 변경 없음. WSL Python3.12.3 / PySCF2.14.0 / ASE3.29.0 / NumPy2.5.4. Windows GUI 시험은 `QUANTUMSTARTER_TK_HIDDEN=1`로 루트·Toplevel 모두 숨김; 화면 배치의 시각 검사는 하지 않았다.

계산·시험 실행은 2026-10-30, 최종 저장 확인·완료 기록은 2026-11-01. 최종 Windows CLI(`main.py --config examples/silicon_cell.json --no-show`)도 exit0/7단계 수렴: `outputs/cell_validation/final_cli/20261030_202319_940271_g2izg_ts/`. version0.9.3, 최초 요청 설정·`final/` 독립 재계산·최종 구조·단계 기록·그림을 확인했다.

**수치 방법.** PySCF FFTDF+동일 uniform XC mesh의 힘·응력과 ASE BFGS/FrechetCellFilter. 원자 및 셀 동시 완화, 외부 압력0. 고정 정수 FFT mesh와 상대 k점. smearing이면 F=E−TS, 없으면 E를 최소화한다. 최대 원자 힘0.02 eV/Å, 응력 잔차0.1 GPa, 단계 E/F 변화1e-5 eV/cell을 각각 검사하고 최종 구조를 다시 SCF 계산한다. 응력 잔차는 isotropic의 |tr(σ)/3|, full의 max|sym(σ)|이다. PySCF 2.14 응력 API는 운동에너지·겹침 적분 내부 유한차분(strain1e-5)을 포함한다. 힘은 해석적 원자 기울기이며 전체 응력이 순수 해석적이라는 이전 문구를 정정했다.

**실제 셀 완화.** Si2, PBE/GTH-SZV-MOLOPT-SR/GTH-PBE, 2×2×2 k점, FFT25³, precision1e-8, SCF1e-9 Ha, 1스레드. 시작 격자 a=5.43 Å의 원시 셀에서 두 번째 원자 x+0.08 Å. full 사례는 추가 대칭 xy strain0.03. 단계 최대40.

| 방식 / 목적함수 | 단계 | 최종 체적 [Å³] | 최대 힘 [eV/Å] | 응력 잔차 [GPa] | 최종 E/F [eV/cell] |
| --- | ---: | ---: | ---: | ---: | ---: |
| isotropic / F (smearing0.1 eV) | 7 | 45.808030 | 1.44886e-4 | 7.48862e-4 | -210.5274642378 |
| full / E (smearing0) | 10 | 45.806122 | 2.66441e-4 | 8.71519e-3 | -210.5274605874 |

최종 독립 SCF와 마지막 단계 목적함수 차이는 두 사례 모두 출력 정밀도에서0. 최종 구조의 모든 길이(좌표와 셀)를 0.995/1.005배로 변형한 독립 단일점에서 목적함수는 각각 isotropic +0.00260741/+0.00252960 eV, full +0.00262046/+0.00251840 eV 상승했다. 이 비교는 샘플링한 체적 방향의 국소 점검이며 모든 변위의 최소점/포논 안정성 증명은 아니다.

원시 결과: `outputs/cell_validation/si_k2/`, `si_full_shear/`의 `validation.json`, `relax/steps/`, `volume_*`. 초기 개발 실행은 버전 필드를 올리기 전이므로 내부 엔진 version이0.9.2이고 최종 재계산을 루트에 저장한 이전 출력 형식이다. 최종 코드는 `final/`에서 재검사한 뒤 루트에 전체 판정을 쓰며 요청 `settings.json`을 보존한다. 최종 저장 방식의 성공·미수렴·취소는 아래 자동 검사와 `final_cli/` 실제 CLI로 확인한다.

**FFT mesh 검사.** 같은 Si 입력 및 isotropic 최종 구조에서 21³→25³→29³. 마지막 두 구간의 모든 차이가 기준(1 meV/atom, 0.01 eV/Å, 0.1 GPa)을 통과했다. 마지막29³를 선택으로 기록하며 원래 입력 mesh를 바꾸지 않는다.

| 구조 / 구간 | ΔF [eV/atom] | 최대 힘 차이 [eV/Å] | 최대 응력 차이 [GPa] |
| --- | ---: | ---: | ---: |
| 입력 21³→25³ | 5.62323e-6 | 3.67854e-5 | 5.45175e-4 |
| 입력 25³→29³ | 1.43959e-6 | 8.26827e-6 | 9.05713e-5 |
| 최종 21³→25³ | 3.35678e-6 | 1.97512e-4 | 1.06727e-3 |
| 최종 25³→29³ | 6.61334e-7 | 5.62373e-5 | 2.77383e-4 |

원시 JSON/CSV와 개별 실행: `outputs/cell_validation/mesh_initial/20261030_201133_825216_psp1okj6/`, `mesh_final/20261030_201518_785695_1jh4k0ut/`. 수렴검사 대상은 고정 구조이며 k점·기저 수렴 및 더 촘촘한 전체 mesh 범위의 보증은 아니다.

**실패와 비교 사례.** 첫 Γ/FFT25³/smearing0 실행은 step3 SCF 미수렴으로 중단됐고 힘을 최적화에 넘기지 않았다(`first_si/`). Γ/smearing0.1은21단계로 수렴했지만 체적 약131.88 Å³, 큰 편차/전단 응력이 남았다(`smeared_si/`); 연구용 Si 평형 자료로 사용하지 않는다. 이 결과와 2³ k점 결과의 큰 차이는 k점 수렴을 FFT mesh 수렴과 분리해야 하는 이유다. 예제의 최소 기저와 낮은 k점 수는 구현 검증용이다.

**자동 검사.** WSL 고유51개, Windows 숨김 GUI6개 통과.

- 핵심27개: `test_periodic_cell`, `test_periodic_convergence`, `test_periodic_gradient`, `test_periodic`. 독립 탄성 모형의 전단/체적/원자 최소점, 자동 mesh 최초1회 결정, F 선택, 단계 한도, 최종 재계산 불일치·취소, 실제 작업자 SCF 실패·단계 사이 취소, CLI 옵션, 힘·응력 에너지 유한차분, 기존 slab/스핀/후보 검사. `outputs/cell_validation/core_regression/tests.log` 및 `tests.json`, 65.286초.
- 밴드/DOS·에너지 분석8개: `test_periodic_analysis`, `test_energy_analysis`. 완화 후 최종 셀이 다른 결과의 에너지 합산 차단 포함. `remaining_tests/tests.log` 및 `tests.json`, 162.337초.
- 공유 실행 회귀17개(위 자동 mesh1개 중복): 분자 입력8, 자체/외부 체크포인트2, 분자 CLI 미수렴1, 수소 저장5, 자동 mesh1. 29.430초, exit0. 물리 허용오차를 변경하지 않았다.
- 최종 GUI6개: 설정·밴드/DOS 필드 보존, FFT 보고 표시, 실제 셀 한도 종료, 작업자 시작 전 취소, 기존 GDF 실계산·그림. `gui_final/tests.log` 및 `tests.json`, 26.300초. 별도 창이 나타나지 않도록 숨김 모드만 사용.
- 첫 회귀 묶음은 Si 밴드 시험에서 요약 없이 exit1로 종료되고 WSL Python 프로세스가 사라졌다. 원인은 확정하지 못했으며 통과로 세지 않았다. 파일 로그/실제 반환 코드를 남기는 runner로 핵심27개와 분석8개를 분리해 재검증해 모두 exit0. 외부 ASE/NumPy deprecation 및 extxyz metadata 경고는 남아 있다.

재현(레이블은 기존 출력과 겹치지 않게 지정):

```powershell
wsl --exec .venv-dft/bin/python scripts/validate_cell.py --label new_isotropic --check-volumes
wsl --exec .venv-dft/bin/python scripts/validate_cell.py --label new_full --mode full --shear .03 --smearing 0 --check-volumes
wsl --exec .venv-dft/bin/python scripts/run_validation_tests.py --label new_core test_periodic_cell test_periodic_convergence test_periodic_gradient test_periodic
```

사용법/제한: [PERIODIC_CELL_KO.md](PERIODIC_CELL_KO.md). slab/GDF 힘·셀 완화, 외부 압력, 제약·최적화 재시작, 포논 안정성·자동 k점/기저 수렴은 미지원. 장기 요청의 IRC/quasi-RRHO/용매/외부 SCF 안정성은 이번 단계 범위 밖이다.

## 0.9.2: bulk FFTDF 힘·응력 제품 API

작업 ID: QS-20261029-NATIVE-DFT 단계 11 (담당 Claude, 2026-10-30; Codex 단계 10 인계 이어받음). 환경 변경 없음. GUI 시험은 숨김 모드만.

- 사전 실측(Si 2원자 한 원자 x 변위, GTH-SZV-MOLOPT-SR, PBE, FFT mesh 25³, UniformGrids, conv_tol 1e-11; scratchpad 프로브): UKS(spin 0) 기울기 vs 에너지 FD 7.4e-9 Ha/bohr, 응력 xx vs 변형 FD 9e-12 Ha/bohr³; 2×2×2 k 메쉬 4.6e-9 / 3.3e-11; Fermi smearing 0.1 eV에서는 기울기·응력이 **자유에너지 F의 미분**과 일치(1e-8 / 1e-11)하고 내부에너지 E 미분과는 1.8e-4 차. → 세 경우 모두 지원, 결과에 미분 대상을 명시.
- `test_periodic_gradient` **4개 통과(61.7초)**: 제품 API(`calculate`)로 기울기 vs FD(허용 1e-7), 응력 xx vs 변형 FD(1e-9), eV/Å·GPa 변환(1 Ha/bohr³ = 29421.0157 GPa), 압력 = −tr/3, 합력 < 1e-3 eV/Å, smearing은 F 미분, GDF와 FFTDF 실행의 반응 분석 혼합 거부. 첫 실행의 입력 시험 1개 실패는 시험 사례가 calculation=gradient를 유지해 다른 규칙이 먼저 걸린 시험 쪽 오류였고 사례를 고쳐 통과.
- 회귀: WSL test_periodic·periodic_analysis·energy_analysis·basis_overrides·dispersion 30개(181.7초), native·molecular·gradient·optimize·frequencies·stability·native_response 52개(101.0초) 통과. Windows(숨김 GUI) 33개 중 30 통과·3 skip(WSL 통과). NEB·제약·Newton 시험은 변경과 무관해 재실행하지 않았다.

한계: bulk 전용(slab은 GDF), FFT mesh 수렴은 자동 확인하지 않음, 셀 최적화 미구현, 시험 구조(변위 Si, SZV 최소 기저)의 응력 값 자체는 물리적 의미가 없고 구현 일치 검사용이다.

## 0.9.1: Cr·Co Newton 수렴, 해석적 응답, 주기계 후속 검토

작업 ID: QS-20261029-NATIVE-DFT 단계 10, ChatGPT Codex, 2026-10-30. 시작 HEAD `8023fc6`. 환경·물리 모델·정수 점유·수렴 허용오차 변경 없음. GUI는 숨김 모드만 실행.

**자체 SCF 및 응답 검증.**
- 기존 Cr(PBE/STO-3G, 2S=6, grid3, 1스레드) 50회: 144.8초, 미수렴, 잔차3.9276e-5. `outputs/newton_validation/baseline_20261030/Cr/report.json`에 이력 보존. 기존 FD 응답/무전처리20회CG/느슨한 에너지 수락 조건의 정체를 재현했다.
- `native_response.py` 해석적 J/K/XC 응답, 궤도 간격 전처리, 적응 CG 정밀도, 반올림 수준의 기울기 개선 검사, 최종 정준화 뒤 전체 기준 재검사. XC 함수 2차 미분만 Libxc에서 가져오며 PySCF SCF/응답 조립은 호출하지 않는다.
- `test_native_response`: RKS/UKS × LDA/VWN·PBE·PBE0·B3LYP 8조건. Fock/궤도기울기 유한차분, 독립 PySCF 응답, H·v 대칭성·선형성, 캐시 및 스트리밍 경로를 검증. 기존 `test_stability`는 에너지 유한차분도 확인한다. core 궤도의 극단 스핀 밀도에서는 FD 비매끄러움(~1e-6)이 관찰돼 작은 일반 궤도 회전의 비정상점에서 검사하며 허용오차는 완화하지 않았다.
- Cr·Co의 수렴 밀도에서 독립 PySCF Newton 재시작도 수렴하고 에너지가1e-9 Ha 이내 일치. 다른 초기 해가 동일하다는 주장이나 전역 바닥상태 보장은 아니다. OH 안장점보다1e-7 Ha 이상 낮은 안정해, 물 DIIS 일치, 정수 점유·전자 수·0/1차원 안정성·강제 미수렴도 검사.

| 원자 | 2S | 반복 | 에너지 [Ha] | 최종 잔차 | 최저 실수 내부 곡률 [Ha] |
| --- | ---: | ---: | ---: | ---: | ---: |
| Cr | 6 | 73 | -1032.8088678851316 | 2.10e-13 | 9.51e-7 |
| Co | 3 | 52 | -1367.1350701366355 | 1.03e-12 | 2.36e-4 |
| Fe | 4 | 53 | -1249.5569908940256 | 7.24e-13 | 4.31e-4 |

공통 PBE/STO-3G/grid3/1스레드, 에너지1e-9 Ha·밀도1e-7·잔차1e-6. Cr·Co 기록 `outputs/newton_validation/analytic_pcg_20261030/`, Fe는 `benchmark_20261030/Fe/`. 전체 SCF 측정117.37/64.68/50.42초는 다른 검증과 동시 실행했으므로 엄격한 성능 비교에 쓰지 않는다.

**같은 수렴 궤도의 응답 속도.** 무작위 방향 고정(seed47), 각7회 중앙값, 한 스레드. 물: 해석0.032245초 vs FD0.120587초(3.740배), 최대 H·v 차3.04e-8. Fe: 해석0.028735초 vs FD0.107455초(3.740배), 차1.62e-8. 최초 캐시 구축 포함 응답은 물0.07681초/Fe0.11395초. 전체 SCF 배속을 의미하지 않는다. 재현: `wsl --exec .venv-dft/bin/python scripts/validate_newton.py --case Fe --max-cycle 100 --label <새이름> --benchmark`.

**slab 밴드.** graphene, a=2.46 Å·셀 높이16 Å, PBE/GTH-SZV-MOLOPT-SR, 2×2×1, smearing0.05 eV, grid2/precision1e-7, SCF1e-10, GMKG16점. 수렴, kz=0, 메쉬 고유값 차3.80e-9 eV, DOS 적분 오차3.31e-7 states, 전자 수 오차1.78e-15. 원시 결과 `outputs/periodic_extension_validation/stage10_20261030/slab/graphene/`, 회귀는 `test_periodic_analysis`. 물성/격자/k메쉬 수렴 검증은 아니다.

**bulk FFTDF 힘·응력 검토.** 변위 Si(두 번째 원자 x+0.08 Å), PBE/GTH-SZV-MOLOPT-SR, Γ, uniform 고정 FFT25³, precision1e-8, SCF1e-11, smearing 없음. 힘의 반대 부호인 에너지 기울기를 원자 변위±0.001 bohr 에너지 FD와 비교: 오차7.39e-9 Ha/bohr. 응력 정의 σ=(1/V)∂E/∂ε, 고정 분율 좌표·±1e-4 변형에서 xx 오차9.27e-12, **비영 yz 전단** 오차1.25e-12 Ha/bohr³(σ_yz=0.0003806094765). 해당 PySCF 응력 구현은 일부 항에 내부 유한차분을 사용하므로 완전 해석적 구현이라고 부르지 않는다. `scripts/validate_periodic_extensions.py fftdf --label <새이름>`, 결과 `outputs/periodic_extension_validation/stage10_shear_20261030/fftdf/report.json`. 이 스크립트는 지원 검토이며 제품 API가 아니다. 일반 k점/UKS/smearing/메쉬 수렴 및 입력·출력·GUI 통합은 후속이다.

회귀: **WSL 고유69개 확인**(첫 관련54개 + 새 수치 경계2개 + 주기계13개; 중복 실행 제외). He 회전 차원0의 나눗셈 경고를 0잔차 조기 종료로 수정한 최종 native/Newton/응답/안정성10개는 모든 경고=오류에서도 통과했다. 최종23개 실행 중15개 통과,8개는 설치된 ASE `atoms.py`의 NumPy2.5 shape 할당 DeprecationWarning 때문에 오류로 종료했다. 해당8개만 수치 RuntimeWarning=오류 조건으로 재실행해 모두 통과(24.802초); 환경/라이브러리 수정은 하지 않았고 외부 deprecation 경고는 남는다. **Windows 숨김 GUI·workflow·수치 경계21개 통과**, skip 없음(최종77.997초, subprocess 실제 returncode=0). 첫 Windows 검사는 OK였으나 PowerShell stderr 리다이렉트의 종료코드1 때문에 실제 프로세스 코드를 수집하여 재검증했다. 원문 로그: `outputs/newton_validation/regression_wsl.log`, `final_focus_wsl.log`, `periodic_followup_wsl.log`, `windows_verified.log`, `windows_exit_status.json`. Python 컴파일·Git 공백 검사·변경 파일 LF 확인 통과.

미완료: bulk FFTDF 제품 통합, IRC, quasi-RRHO/용매, 삼중항·복소 안정성, 자동 바닥상태 탐색. HEA는 현재 연구 목표가 아니므로 과거 예시 검증을 새 우선순위로 채택하지 않았다.

## 0.9.0: 자체 2차 SCF, 주기계 밴드·DOS

작업 ID: QS-20261029-NATIVE-DFT 단계 9 (담당 Claude, 2026-10-30). 환경 변경 없음. GUI 시험은 숨김 모드만.

**자체 Newton SCF (engines/native_newton.py).**
- `test_newton` 2개 통과(245.6초, 백그라운드 Cr 계산과 동시 실행): 물(PBE/STO-3G/grid2) Newton = DIIS 에너지 1e-8 이내, 최종 잔차 < 기준, aufbau 위반 없음; OH(grid1 엄격 수렴)에서 DIIS가 수렴한 불안정 해(최저 고유값 < −1e-4)에서 출발해 탈출 기록을 남기고 더 낮은(>1e-7 Ha) 안정 최소점에서 종료, 전자 수 9.
- 프로브: 물 def2-SVP 9회(36.8초)·Fe STO-3G 57회(183초), 둘 다 DIIS와 같은 에너지. 탈출 기능 이전 Cr은 안장점(최저 고유값 −0.42, 0.13 Ha 높음)에서 '수렴' 판정 → 수렴 시 고유값 검사와 탈출을 추가. 이후 Cr 100회·300회 모두 미수렴(300회: 667초, 기울기 4~7e-5 진동, E −1032.8088625, PySCF DIIS 수렴해 −1032.8088486보다 1.4e-5 Ha 낮음). Co는 시간 제한으로 결과 없음.
- 개발 중 버그: pyscf_engine이 scf_solver=newton이면 자체 엔진에도 PySCF `mf.newton()` 래퍼와 `env['imacro']`를 적용 → PySCF 엔진에만 적용하도록 수정. PySCF Newton 기준 비교 프로브는 내 래퍼 사용 오류(물 −75.795 Ha)로 무효 처리하고 PySCF DIIS 값을 기준으로 썼다.

**주기계 DOS·밴드 (engines/periodic_analysis.py).**
- PySCF 2.14 실측: GDF 주기계 기울기 미지원(`'GDF' object has no attribute 'get_j_e1'`), 응력은 FFTDF 전용(AssertionError) → 힘·응력은 구현하지 않음. `get_bands`는 GDF에서 동작.
- Si(examples/silicon_bulk.json, GTH-DZVP-MOLOPT-SR, PBE) 2×2×2, smearing 0, conv_tol 1e-10, 경로 auto(GXWKGLUWLK,UX) 40점: 253초. DOS 적분 51.999999(기대 52), 가중 점유 합 8 = 원자가 전자, 메쉬 HOMO/LUMO 갭 0.616 eV, 경로 간접 갭 0.488 eV(VBM Γ, CBM Γ–X 위), 직접 갭(Γ) 2.51 eV, 메쉬 k점 밴드값과 SCF 고유값 차 1.1 meV(conv_tol 1e-10에서도 동일 → get_bands 경로의 고유 차이로 판단). 기저·k점 수렴 미확인이므로 정량 갭 값으로 쓰지 않는다.
- `test_periodic_analysis` 3개 통과: 합성 고유값에서 DOS 적분·점유 합·HOMO/LUMO(RKS·UKS·smearing 금속), Γ점 Si 실계산의 DOS 적분·점유 합·경로 라벨·메쉬 일관성 < 5 meV·직접 ≥ 간접 갭.
- 회귀: WSL 89개 통과(68개 121.8초 + 21개 230.6초). Windows(숨김 GUI) 31개 중 28 통과·3 skip(WSL 통과).

한계: Newton은 느리고 Cr·Co 미해결; 주기계 힘·응력 없음; slab 밴드 경로 실계산 미검증; 밴드는 비자기일관 고유값이며 SOC 없음.

## 0.8.9: CI-NEB

작업 ID: QS-20261029-NATIVE-DFT 단계 8 (담당 Claude, 2026-10-30). 환경 변경 없음. GUI 시험은 숨김 모드만.

- 확인한 geomeTRIC 1.1.1 동작: NEB는 입력 프레임을 이미지로 쓰고(보간 없음, `images`는 균등 추출), 양 끝은 고정(locks), 수렴 maxg/avgg 단위는 eV/Å(기본 0.05/0.025), `OptimizeChain`은 최대 반복에 도달해도 반환하므로 반환 후 기준 충족을 직접 판정한다.
- `test_neb` **3개 통과(25.2초)**: NH₃(PBE/STO-3G/grid3, 자체 엔진) tight 최적화 피라미드 ↔ z 거울상, 7 이미지. 개발 실행: 5회 수렴, 기울기 계산 27회(내부 5 × 5 + 양 끝 2), 최대/평균 힘 0.0258/0.0132 eV/Å, 상대 에너지 0, 17.969, 45.833, **54.986**, 45.833, 17.969, 0 kJ/mol(거울 대칭), 최고 이미지 평면(2.6e-8 Å), 그 구조에서 transition_state 정밀화 장벽 54.985 kJ/mol(차 0.001 kJ/mol; 시험 허용 0.1), 양 끝 에너지 차 1.8e-9 Ha. 기록 파일 수 검사, geomeTRIC 부재 시 명확한 오류.
- 개발 중 문제: Windows Python `Path.write_text`가 수정한 소스 7개를 CRLF로 바꿈(첫 CRLF 검사 명령 `grep -c $'\r$'`은 Git Bash에서 모든 줄을 세는 오탐이었고 바이트 검사로 확인). 커밋 전 LF로 복원, 이후 스크립트 수정은 bytes로 쓴다.
- 회귀: WSL 84개 통과(136.5초). Windows(숨김 GUI) 36개 중 32 통과·4 skip(WSL 통과).

한계: 초기 경로는 선형 보간뿐(IDPP 등 없음), 자동 정렬 없음, IRC 없음, 하나의 경로만 탐색(다른 반응 경로 가능성 배제 못 함).

## 0.8.8: 제약 최적화·이완 스캔·전이상태 탐색

작업 ID: QS-20261029-NATIVE-DFT 단계 7-C (담당 Claude). 환경 변경 없음. GUI 시험은 숨김 모드만.

- `wsl ... -m unittest test_constraints -v`: **4개 통과(30.0초)**. 전부 자체 엔진 PBE/STO-3G/grid3.
  - 물 `$set angle 2 1 3 100`: 최종 각도 100.00000000000003°(시험 1e-6), 두 O–H 1.035081 Å 대칭, 제약 하에서 힘 기준 판정은 기록하지 않음(None).
  - 물 `$freeze xyz 1`: O 원자 이동 3.4e-5 Å(시험 < 1e-3 Å). geomeTRIC은 제약을 수렴 허용오차 안에서 강제하며 정확히 0은 아니다.
  - 잘못된 제약(`angle 2 1 100`) → 실행 전 parse_constraints 단계에서 ValueError.
  - H₂ 이완 스캔 0.6–1.0 Å 5점: 거리 오차 < 1e-6 Å, 같은 거리 단일점 에너지와 2e-10 Ha 이내(이원자는 다른 자유도가 없으므로 같아야 함), 최저점 0.7 Å, scan.csv·상단 density.cube 생성.
  - NH₃ 반전 전이상태(거의 평면 시작): 29회 기울기 계산(초기 Hessian 24 = 4원자×3×2, geomeTRIC dirname의 hessian 표시로 분류), 평면 구조(N과 H 평균 z 차 1.1e-4 Å), **진동수: 허수 정확히 1개 −880.6 cm⁻¹, 분자면 수직 우산 모드**, 판정 transition_state_candidate.
- 개발 중 오류: 사용자 제약 문자열이 최적화 내부의 단계별 기울기 설정(calculation=gradient)으로 넘어가 새 검증 규칙에 걸림 → 내부 단계·최종 재계산 설정에서 제약을 비움.
- 외부 라이브러리 경고: geomeTRIC 1.1.1 `optimize.py`가 제약 파일을 닫지 않아 ResourceWarning(결과 무관).
- 회귀: WSL 81개 통과(120.0초). Windows(숨김 GUI) test_molecular_gui·test_workflow·test_solver·test_constraints 30개 중 27 통과·3 skip(WSL 통과).

한계: NEB·IRC 없음(geomeTRIC 1.1.1에 `ElasticBand`/`OptimizeChain`/`NEBParams` 존재 확인 — 다음 단계 후보), 스캔은 1차원만, 전이상태는 시작 구조 근처의 1차 안장점만 찾으며 원하는 반응 경로인지 사용자가 확인해야 한다.

## 0.8.7: 이상기체 자유에너지 (RRHO)

작업 ID: QS-20261029-NATIVE-DFT 단계 7-A (담당 Claude). 환경 변경 없음. GUI 시험은 숨김 모드만.

- 확인한 사실: PySCF `hessian.thermo.thermo`는 병진·회전 항에 동위원소 평균 질량을 쓴다. 우리 결과는 모든 항에서 최빈 동위원소 질량(Gaussian 관례). 비교 시에는 같은 평균 질량을 우리 함수에 넣었다.
- 같은 Hessian·질량·진동수에서 PySCF thermo와 항목별 비교(물, CO, NH₃, 삼중항 O₂; 개발 확인): S_trans·S_rot·S_elec 상대 차 3.3~4.2e-7(물리상수 CODATA 판 차이와 같은 크기의 일정한 차), ZPE 1e-10, H·G 상대 1e-10 이하(≈1e-8 Ha). 선형/비선형 판정과 대칭수(물 2, CO 1, O₂ 2) 일치.
- `test_frequencies` **8개 통과(26.3초, RuntimeWarning을 오류로 둔 실행)**: 위 비교(엔트로피 상대 1e-5, H·G 1e-7 Ha), 반응 분석에서 G 선택 시 값 = 수동 합, 원자의 rotor = atom, 온도가 다른 항(298.15 vs 500 K) 거부, 같은 실행의 전자 에너지 반응은 온도와 무관하게 허용.
- 개발 중 경고: 원자(관성 모멘트 0)에서 쓰이지 않는 회전 온도 계산의 0 나눗셈 경고 → 회전체일 때만 계산하도록 수정.
- 회귀: WSL 77개 통과(test_frequencies 8 + 나머지 69). Windows(숨김 GUI 포함) 37개 중 31 통과·6 skip(PySCF 필요, WSL 통과).

한계: RRHO만(quasi-RRHO·비조화·내부 회전자 없음), 들뜬 전자 상태 무시, 용매·표준상태(1 M) 변환 없음, 대칭수 자동 판정은 PySCF 점군 판정의 허용오차에 의존(필요하면 `symmetry_number`로 지정).

## 0.8.6: 자체 SCF 품질 — 초기 추정·혼합·안정성 분석

작업 ID: QS-20261029-NATIVE-DFT 단계 7-B (담당 Claude). 환경 변경 없음. **사용자 요청에 따라 이번부터 GUI 시험은 숨김 모드(`QUANTUMSTARTER_TK_HIDDEN=1`)로만 실행**했다.

- 반복 수 진단(PBE·B3LYP, def2-SVP, grid3): 기본 혼합 0.15가 쉬운 계에 3~4회, core 추정이 2~3회 추가. 혼합 0은 Fe 원자 미수렴, 항상 혼합은 쉬운 계 손해 → 잔차 기준 혼합. 기준 1e-2/1e-3/1e-4/0에서 물 9/9/9/12, CH4 7/7/8/11, Fe 173/141/미수렴/미수렴회 → 1e-3 채택.
- 첫 구현은 minao 밀도를 SCF 시작 밀도로 직접 썼고 test_native의 취소→재시작 시험이 "체크포인트 스핀별 전자 수 불일치"로 실패(투영 밀도의 전자 수가 N이 아님 + 첫 반복 혼합). 첫 Fock 생성에만 쓰도록 수정 후 통과.
- 결과(기본 설정): 물 14→9, CH4 14→7, NH2 14→10, O2 10→8회, 에너지 불변(≤1e-11 Ha). Fe(STO-3G, 2S=4) 미수렴→수렴(test_native 회귀 시험).
- Cr(2S=6)·Co 원자(STO-3G/def2-SVP): 모든 조합 미수렴. 페르미 스미어링 σ = 0.01/0.03 Ha 상태는 13~25회 수렴하나 정수 점유로 되돌리면 진동. 진동 중 에너지가 PySCF 수렴해보다 1e-5 Ha 낮음 → PySCF 해도 최저가 아닐 수 있음. 미해결.
- 안정성 분석 `test_stability` **3개 통과(26.3초)**: 비정지점(core 궤도)에서 기울기 공식 = 에너지 중심차분(상대 1e-6), 수렴 해에서 x·Hx = 에너지 2차 차분(상대 1e-4; 개발 확인 물 15.03290/15.03290, OH 2.404461/2.404468), 물(PBE0) 안정·OH(UKS PBE) 불안정 판정이 PySCF `stability()`와 일치, follow가 OH에서 9.3e-7 Ha 낮은 해로 이동하고 최종 최저 고유값 +2.2e-4(안정). 이 해는 PySCF가 찾은 해보다도 낮다.
- 회귀: WSL 75개 통과(94.0초). Windows(숨김 GUI) test_molecular_gui·test_workflow·test_solver·test_stability 29개 중 27 통과·2 skip(PySCF 필요, WSL 통과). 수소·주기계 GUI 시험은 숨김 모드가 없어 실행하지 않았다.

## 0.8.5: 조화 진동수·ZPE

작업 ID: QS-20261029-NATIVE-DFT 단계 6 (담당 Claude). 환경 변경 없음.

- 질량: PySCF `atom_mass_list(isotope_avg=False)`는 질량수 정수(16, 1, 1)를 반환함을 확인했다(H 진동수 ~0.4% 오차 원인). 사용하지 않고 `pyscf.data.elements.COMMON_ISOTOPE_MASSES`(1.007825, 15.994915 …)를 쓴다.
- `wsl ... -m unittest test_frequencies -v`: **6개 통과(25.8초)**.
  - 투영·질량가중·단위 독립 확인: 같은 Hessian에 대해 우리 `harmonic_analysis`와 PySCF `thermo.harmonic_analysis`가 물·CO(선형)·NH₃에서 1e-4 cm⁻¹ 이내 일치.
  - 물 PBE/STO-3G/grid3(tight 최적화 구조): PySCF 해석적 Hessian 진동수와 2 cm⁻¹ 이내, 허수 없음, minimum, 변위 폴더 19개, ZPE = ½Σν.
  - H₂: 1개 모드(선형), 에너지 곡률 d²E/dr²에서 얻은 진동수와 1 cm⁻¹ 이내.
  - 평면 NH₃: |ν| > 20 cm⁻¹ 허수 진동 ≥1, 그 모드는 분자면 수직(z) 변위, 경고 포함. 변위 계산 SCF 실패 시 중단.
  - 열역학 공식: T = 1 K에서 유한값(경고를 오류로 둔 실행), 고온에서 모드당 Cv → k_B.
- 개발 비교(PBE/def2-SVP/grid3, tight 최적화 구조): 자체 1608.59/3691.20/3790.51 vs PySCF 해석적 Hessian 1608.64/3690.72/3790.50 cm⁻¹; 우리 Hessian을 PySCF harmonic_analysis에 넣은 값은 소수 3자리 동일. H₂ 4277.63 vs d²E/dr² 4277.72 cm⁻¹.
- 실제 CLI(Windows→WSL) 물 PBE-D4/def2-SVP/grid3, 0.8.4 CLI 최적화 구조 입력, `outputs/frequency_validation/cli/`: 자체 1608.48/3690.73/3790.01, PySCF 1608.43/3690.50/3789.62 cm⁻¹, 둘 다 minimum, ZPE 54.366/54.362 kJ/mol, 병진·회전 잔여 최대 16.6/17.1 cm⁻¹, Hessian 비대칭 4.3e-5/6.2e-5 Ha/bohr². 실행 당시 version 0.8.4.
- 오류·조치: (1) T = 1 K에서 expm1 오버플로로 Cv가 NaN이 될 수 있어(JSON 저장 실패 위험) e^(−x) 기반 안정식으로 수정. (2) Windows 전체 실행에서 **기존 수소 GUI 시험이 간헐 실패**(3회 중 2회): 실행 폴더 이름이 초 단위 시각 + 무작위 문자열이라 같은 1초 안의 두 실행 순서가 무작위로 정렬됨(기록 목록의 실제 버그). 수소·분자 실행 폴더 이름에 마이크로초를 넣어 수정, 실패 조합 5회 반복 통과. 처음에 진동수 GUI 시험과 상관이 있어 보였던 것은 1회씩 시도한 우연으로 판단.
- 회귀: WSL 71개 통과(82.9초). Windows 63개 중 44 통과·19 skip(18개는 PySCF/dftd4/geomeTRIC 필요로 WSL 통과, 기존 주기계 GUI 1개).

한계: 조화 근사만(비조화·보정 인자 없음), 병진·회전·전자 열역학과 G 없음, 유한차분 Hessian은 6N회 계산이 필요해 큰 분자에서 느리다(해석적 Hessian 없음).

## 0.8.4: 자체 DIIS 정체 수정과 구조 최적화

작업 ID: QS-20261029-NATIVE-DFT 단계 5 (담당 Claude). WSL `.venv-dft`에 geomeTRIC 1.1.1(+networkx 3.7) 추가 설치.

**DIIS 수정 (커밋 de7d245).** OH(UKS PBE/STO-3G, 격자 1, 엄격 수렴)가 300회 미수렴. 진단: β 점유공간 연속 겹침 0.9999999999·HOMO-LUMO 간격 0.027 Ha 일정 → 점유 교대 아님.
Pulay 계수 상한 30 때문에 300회 중 289회 이력이 버려졌다(밀도 혼합 0, 레벨 시프트 0.1/0.2, DIIS 16으로는 불변). 의사역행렬 풀이로 바꾼 뒤 OH 25회, 계수 상한 1e3도 22회 수렴.
H2O(def2-SVP)/NH2/O2/선형 OH의 반복 수와 에너지는 변화 없음(1e-13 Ha). 회귀 시험은 수정 전 코드에서 실패함을 확인했다.
수렴한 OH 해는 PySCF 해보다 +8.27e-6 Ha이며 PySCF `stability()`가 내부 불안정으로 판정(불안정 방향 추종 시 −7.98e-6 Ha). PySCF 해도 불안정 판정 → 다중 UKS 해. Fe 원자(STO-3G, 2S=4)는 수정 전후 모두 미수렴.

**구조 최적화.**
- `wsl ... -m unittest test_optimize -v`: **4개 통과, 19.6초**. 물 PBE/STO-3G/grid3/tight: 자체 엔진+우리 geomeTRIC 연결 vs **PySCF 자체 SCF·기울기·geomeTRIC 드라이버(pyscf.geomopt)** 구조 차 < 2e-4 Å/°, 최종 최대 기울기 < 기준, 기록 파일·단계 폴더·이전 구조 밀도 초기 추정 확인, 최대 단계 1 → not_converged, 단계 SCF 실패 시 중단, 중단 요청 반영, geomeTRIC 부재 시 명확한 오류, 작업자 프로세스 경유 실행.
- 개발 중 비교(tight, PBE/def2-SVP, grid 3): H2 자체 6단계 r = 0.768334 Å, H2O 자체·PySCF 각 5단계 r(OH) = 0.974677 Å, ∠HOH = 102.044013° — 독립 드라이버와 소수 6자리 일치. 최종 재계산 에너지 차 ≤ 2.4e-10 Ha.
- 실제 CLI(Windows→WSL) 왜곡 물 `examples/water_distorted.xyz`, PBE-D4/def2-SVP/grid3/default, `outputs/optimization_validation/cli/`: 두 엔진 모두 4단계 수렴, E = −76.2727809465 Ha, ΔE_opt = −0.048490 eV, r(OH) = 0.974723 Å, ∠HOH = 102.0410°, 최대 힘 5.6e-4 eV/Å. 단계당 SCF 자체 14/9/10/8회, PySCF 7/6/6/5회. 실행 당시 version 필드 0.8.3.
- 회귀: WSL 65개 통과(75.7초; 첫 실행 1 실패 = test_gradient가 'optimize'를 잘못된 값 예시로 쓴 시험 쪽 문제, 'hessian'으로 교체). Windows 실계산 GUI 포함 56개 중 41 통과·15 skip(PySCF/dftd4/geomeTRIC 필요 14개는 WSL 통과, 기존 주기계 GUI 1개).

한계: 국소 최소점 후보이며 진동수(허수 진동) 확인 없음. 제약·스캔·전이상태 없음. 최적화 체크포인트 재시작 없음. 열린 껍질에서는 단계마다 같은 SCF 해에 머무는지 보장하지 않는다(이전 구조 밀도로 시작해 연속성을 높일 뿐).

## 0.8.3: 해석적 원자 힘

작업 ID: QS-20261029-NATIVE-DFT 단계 4 (담당 Claude). 환경 변경 없음.

- `wsl --exec .venv-dft/bin/python -m unittest test_gradient -v`: **6개 통과, 11초**.
- `wsl --exec .venv-dft/bin/python scripts/validate_gradient.py` → `outputs/gradient_validation/20261029_211817/report.json`:
  - **고정 격자 중심 유한차분**(h = 1e-4 bohr, 격자를 기준 구조에 고정, 모든 원자·축, STO-3G, 격자 1, 엄격 수렴): H₂O LDA 1.17e-9, PBE 3.26e-9, PBE0 2.46e-9, NH₂(UKS) PBE 1.90e-9, B3LYP 1.39e-9 Ha/bohr(최대 |기울기| ~0.1). 고정 격자에서는 해석식이 근사 없는 미분이어야 하므로 이것이 구현 정확성의 주 검증이다.
  - **격자 응답 크기**(격자가 원자와 함께 움직이는 유한차분, h = 1e-3, H₂O PBE): 격자 수준 1/2/3/4에서 최대 차 6.04e-4 / 4.40e-5 / 1.08e-5 / 1.13e-6 Ha/bohr, 합력 8.2e-4 / 3.6e-5 / 1.4e-5 / 3.3e-7. 격자 응답 항을 넣지 않은 근사의 크기이며 격자 수준을 올리면 단조 감소했다.
  - **PySCF 기울기와 비교**(같은 격자, def2-SVP, 격자 3): PBE 8.8e-13, B3LYP 5.7e-13 Ha/bohr.
  - PySCF 경로의 Newton SCF·스칼라 X2C(AgH, Dyall-v2z) 기울기 지원 확인.
- D4 기울기: 유한차분과 1e-9 이내(test_gradient), 두 엔진 결과의 분산 성분이 dftd4 직접 계산과 1e-12 일치.
- AuH(Au def2-ECP): PySCF 경로 기울기 수렴·대칭 확인. SCF 미수렴 시 기울기 없음+경고 확인.
- 실제 CLI(Windows→WSL) 물 PBE-D4/def2-SVP/grid3 `--calculation gradient`, `outputs/gradient_validation/cli/`: 최대 힘 PySCF 1.168150, 자체 1.168146 eV/Å,
  기울기 성분 차 최대 ~1.1e-7 Ha/bohr — 기본 수렴 기준(자체 밀도 1e-7/잔차 1e-6)의 차이이며 엄격 기준에서는 위 1e-12 수준. 합력 8.0e-6 Ha/bohr.
  이 두 실행의 `version` 필드는 0.8.2이다(버전 문자열 변경 직전 실행, 코드는 같음).
- 회귀: WSL 60개 통과(56.7초). Windows 실계산 GUI 포함 50개 중 38 통과·12 skip(PySCF/dftd4 필요 11개는 WSL 통과, 기존 주기계 실계산 GUI 1개). 새 실계산 GUI 기울기 시험 1개 통과(7.8초).

발견한 문제(미해결): OH(UKS, PBE, STO-3G, 격자 1, 밀도 1e-10/잔차 1e-9)에서 자체 SCF가 300회 동안 ΔD ≈ 2e-5에서 진동하며 수렴하지 않았다.
거의 축퇴한 π 궤도의 점유 교대로 추정하며 기울기 구현과 무관하다. 시험은 비축퇴 라디칼 NH₂로 대체했고, 원인 확정과 수정(레벨 시프트·점유 고정 등)은 남은 작업이다.

## 0.8.2: DFT-D4 분산 보정

작업 ID: QS-20261029-NATIVE-DFT 단계 3 (담당 Claude). WSL `.venv-dft`에 dftd4 4.3.1, tomli 2.5.0 추가 설치.
모델: DFT-D4(EEQ 전하 의존 C6, BJ 2체 + ATM 3체), 매개변수는 dftd4가 배포하는 PBE/PBE0/B3LYP 값을 그대로 사용하고 결과에 기록.

- `wsl --exec .venv-dft/bin/python -m unittest test_dispersion -v`: **7개 통과, 9.2초**.
  - 라이브러리 계약: 118원소를 실제로 넣어 지원 집합이 코드의 D4_ELEMENTS(Z 1~103, 112~118)와 같음, LDA 매개변수 없음, 명시 매개변수 = 라이브러리 이름 조회(1e-14).
  - 독립 수식 재현: Ar₂·ArKr에서 2체 D4-BJ 식 E = −C6[s6/(R⁶+f⁶) + s8·rr/(R⁸+f⁸)]을 직접 계산해 R = 3, 10, 50, 200 bohr에서 상대 1e-10 이내 일치. 200 bohr는 라이브러리 기본 cutoff(60 bohr, 61 bohr부터 0으로 실측)를 넘으며 이번 구현이 절단 없이 합산함을 확인. C6 자체(EEQ·CN 의존)는 라이브러리 값을 사용하므로 이 시험의 독립 검증 범위 밖이다.
  - 회전·병진 불변(1e-13), 총전하 전달(EEQ 전하 합 = +1, 에너지 변화), Å/bohr 입력 동일(1e-12).
  - SCF 통합: PySCF/자체 엔진 모두 총에너지 = SCF + 분산(1e-12), 성분 합 일치, 분산을 켜도 SCF 에너지 불변(1e-9), dftd4 패키지 자체 PySCF 어댑터와 1e-12 일치, 두 엔진 총에너지 2e-7 이내.
  - 라이브러리 부재 시 명확한 오류, `dispersion=none`은 정상 수렴.
- 전체 WSL 회귀 `test_native test_molecular test_periodic test_basis_overrides test_energy_analysis test_dispersion`: **54개 통과, 52.4초**.
  첫 실행에서 test_basis_overrides의 "이전 서명 공식" 시험 1개가 실패했다. 시험이 과거 설정을 재현할 때 새 필드 `dispersion`을 빼지 않은 시험 쪽 누락이었고(실제 서명 함수는 dispersion 제외), 시험을 고쳐 재실행 통과.
- Windows `QUANTUMSTARTER_GUI_DFT_TEST=1` `test_molecular_gui test_workflow test_gui test_solver test_periodic_gui test_basis_overrides test_energy_analysis test_dispersion`: **44개 중 37 통과, 7 skip, 61.6초**.
  skip은 PySCF/dftd4가 필요한 6개(WSL에서 통과)와 별도 환경변수가 필요한 기존 주기계 실계산 GUI 1개.
- 실제 CLI(Windows→WSL) 물 PBE/def2-SVP/grid3, `outputs/dispersion_validation/`:
  PySCF+D4 −76.272209814055 Ha(SCF −76.272014472667, D4 −1.953413873068e-4), 자체+D4 −76.272209814057 Ha(차 2e-12),
  D4 없는 PySCF −76.272014472667 Ha(= D4 계산의 SCF 에너지). EEQ 전하 O −0.598, H +0.299.
  이 세 실행의 `version` 필드는 0.8.1이다. 버전 문자열을 0.8.2로 바꾸기 직전에 실행했으며 그 외 코드는 같다.

한계: D4 에너지 정확도(벤치마크 대비)를 이 프로젝트에서 재검증하지 않았다. 매개변수는 해당 범함수 구현으로 피팅된 값이며 PySCF B3LYP(Libxc 402)와 피팅에 쓰인 B3LYP의 VWN 변형이 다를 수 있다.
분산 힘, 주기계 D4, D3 비교는 미구현이다.

## 0.8.1: 원소별 기저 선택

작업 ID: QS-20261029-NATIVE-DFT 단계 2 (담당 Claude). 환경은 0.8과 같다(WSL PySCF 2.14.0, Windows Python 3.13/Tk).

- `wsl --exec .venv-dft/bin/python -m unittest test_basis_overrides test_energy_analysis -v`: **10개 통과, 3.5초**.
  입력 정규화·거부, 원소별 ECP/상대론/격자 판정, 체크포인트 서명(지정이 없으면 이전 공식과 같은 해시), 혼합 기저 실계산, 반응 분석의 원소별 기저 대조.
- 혼합 기저 물(O=def2-SVP, H=STO-3G, LDA/VWN, grid2, conv_tol 1e-10): PySCF 엔진 결과가 PySCF 직접 dict 기저 구성과 1e-8 Ha 이내,
  자체 엔진과 2e-7 Ha 이내, AO 16개(균일 STO-3G 7개). 같은 XC·격자에서 혼합 기저 에너지가 균일 STO-3G보다 낮다(변분 자유도 증가와 일치).
- AuH(Au=def2-SVP+ECP, H=STO-3G, LDA): Au 코어 60, H 0, 명시 전자 20, `Tr(DS)`=20으로 수렴.
- `wsl ... -m unittest test_native test_molecular test_periodic`: **기존 37개 통과, 49.6초** (변경 후 회귀).
- `QUANTUMSTARTER_GUI_DFT_TEST=1` Windows `test_molecular_gui test_workflow test_basis_overrides test_energy_analysis`: **26개 중 24 통과, 2 skip(Windows에 PySCF 없음 — 같은 2개는 위 WSL 실행에서 통과), 54.8초**. 실제 WSL 계산 GUI 5개와 새 원소별 기저 입력 GUI 1개 포함.
- 실제 CLI(Windows→WSL) 혼합 기저 물, LDA/VWN, grid3: PySCF 엔진 -75.770264060981 Ha(6회), 자체 엔진 -75.770264061012 Ha(16회), 차이 3.1e-11 Ha,
  `Tr(DS)`=10. 결과: `outputs/basis_override_validation/20261029_202537_953szltv`, `.../20261029_202543_ocuvwdck`.

한계: 혼합 기저의 총에너지는 다른 기저 조합과 직접 비교할 수 없다. 기저 혼합이 특정 물성에 충분한지는 기저 수렴 비교로 별도 확인해야 하며,
이 검증은 구현 일치(같은 입력에서 같은 해밀토니안)에 대한 것이지 정확도 검증이 아니다. 주기계 엔진은 원소별 기저를 지원하지 않는다.

## 0.8: 자체 분자 SCF와 연구용 결과 분석

작업 ID: QS-20261029-NATIVE-DFT (첫 단계 검증 완료, 첨부 후속 확장은 인계 유지).
환경: WSL PySCF 2.14.0, NumPy/SciPy/Libxc; Windows Python 3.13/Tk. 라이브러리 경계와 계산식은
[자체 엔진 문서](NATIVE_DFT_KO.md), 미완료 항목은 [요구사항 추적](REQUIREMENTS_STATUS_KO.md)에 있다.

- `wsl --exec .venv-dft/bin/python -m unittest test_native test_energy_analysis test_molecular test_periodic -v`: **39개 통과, 66.843초**. 최종 canonical 밀도에도 에너지·밀도·잔차 기준을 적용한 코드에서 실행했다. 원본 출력: `outputs/native_validation/final_physics_tests.log`.
- `$env:QUANTUMSTARTER_GUI_DFT_TEST='1'` 후 Windows `test_molecular_gui -v`: **10개 통과, 64.352초, skip 없음**. 실제 자체 계산 중단→체크포인트 재시작, SCF 그래프 저장, 외부 ECP/X2C/Co Newton, 미수렴·실패 복구·ZIP 내보내기를 포함한다. 원본 출력: `outputs/native_validation/final_gui_tests.log`.
- Windows `test_workflow -v`: **5개 통과, 0.104초**. 기존 수소 검증기의 설정·저장·그림 경로를 보존했다.
- 자체 9조건(He/H2/H2O/H/OH, LDA/PBE/PBE0/B3LYP, STO-3G/def2-SVP)의 계산과 기준 비교. 외부 RKS/UKS kernel 호출을 강제로 실패시키는 상태에서 자체 계산이 통과했다. 이는 외부 적분/XC 의존성까지 제거했다는 뜻은 아니다.
- RKS/UKS×4 XC의 Fock를 에너지 유한차분 미분과 비교(허용오차2e-7), 전자·스핀 수, MO 정규화, 밀도 재구성, 성분 합, 최종 에너지·밀도·잔차 및 이력 일치, 초기 추정/격자, 강제 미수렴, 호환성 오류를 검사했다.
- 영구 5종 비교: `outputs/native_validation/20261029_195245_027252/comparison.json`. 같은 STO-3G/grid2/conv_tol1e-10에서 He/H2/H2O 독립 PySCF 대비 절대차는 각각4.44e-16/4.44e-16/7.11e-14 Ha. OH는 독립 SCF 간2.40551e-6 Ha 차이, 자체 밀도로 외부를 재시작하면2.84e-14 Ha 차이. 서로 다른 비교를 혼동하지 않는다.
- 실제 CLI 자체 물 예제: `outputs/molecular/20261029_195734_a6kkyr2a`, LDA/STO-3G/grid3, -74.7321051890996 Ha, 13회 수렴. 자체 지원은 H~Kr/비상대론 전체전자/점핵이며 이 예제를 전체 원소의 연구 정확도 근거로 쓰지 않는다.
- 실제 분석 H2−2H: -6.799190445 eV, `outputs/analysis/20261029_195626_hr778v3r/analysis.json`. 소스 SHA256·조건·수지를 확인했다. LDA 최소 기저의 합산 검증이며 실험 결합에너지 정확도 검증이 아니다.

발견·수정한 오류와 한계:

1. UKS 정확교환 합산의 NumPy ellipsis 축 오류를 명시적 스핀축 합산으로 수정했다.
2. GUI가 SCF 메시지마다 그래프를 다시 그려 중단 이벤트가 늦어지던 문제를30ms 큐 처리 예산과 배치당1회 갱신으로 수정했다. 중단 시험은 첫 체크포인트 파일을 관찰해 실제 버튼을 누른다.
3. OH의 독립 SCF 에너지 차이를 조사했다. 같은 밀도에서 에너지 차이4.26e-14 Ha/Fock 최대차3.55e-15 Ha이며 그 밀도에서 외부도 같은 에너지로 수렴했다. 가능한 퇴화 궤도·각도 격자·해 선택의 영향을 기록했지만 원인을 확정하거나 전역 바닥상태를 보증하지 않는다.
4. 중단은 안전한 계산 경계에서 반영되며 외부 적분 한 단계 내부에서 즉시 멈추지 않는다. 체크포인트는 밀도를 복원하며 DIIS 이력을 복원하지 않는다.
5. 반응 분석은 수렴·원소/전하 수지·동일 엔진/계산 조건을 검사한다. 서로 다른 주기 셀·k점의 기준 에너지는 현재 거부한다. 진동·열역학·촉매 활성 계산을 수행한 것이 아니다.

아래 0.7 및 이전 기록은 당시 검증 내역으로 보존한다.

## 0.7: 모든 원소의 분자 입력·상대론 경로

작업 ID: QS-20261028-118. PySCF 2.14.0 / WSL Python 3.12.3 / Windows GUI Python 3.13.

- `wsl --exec .venv-dft/bin/python -m unittest test_molecular test_periodic -v`: **30개 모두 통과**, 최종 heavy universal 격자 적용 후 39.159초. 물리적 정확도 성적이 아니라 입력/저장/계산 불변량 및 회귀 검사다.
- 9기저×118원소, 1,062개 조합의 실제 데이터와 프런트 지원표가 일치: 지원554 / 미지원508.
- 두 Dyall 기저의 118원소 각각에 대해 Gaussian 핵·전자 수·AO 수·유한한 겹침 행렬을 검사했다. universal 격자는 118개 원소 모두 생성했다. SCF를 118종 전부 실행한 것은 아니다.
- 모든 지원 def2-SVP/TZVP ECP의 실제 코어 수를 검사했다. 28/46/60개가 원소별 설정과 일치한다.
- X2C의 RKS He, UKS H, DIIS/Newton, 저장·MO 직교성·밀도행렬 재구성·Tr(DS)·스핀, 강제 미수렴 구분, CUBE/그림의 X2C 표기를 검사했다.
- Windows 입력 GUI 기존3개 통과, 원소 지원표/Og 설정 저장1개 통과. `QUANTUMSTARTER_GUI_DFT_TEST=1`로 새 X2C 실계산 GUI1개도 통과(11.339초): Windows→WSL→결과/그림/설정 파일. 이전 세 가지 선택적 실계산 GUI 테스트는 이번 변경에서 재실행하지 않았다.
- `test_workflow` 5개 통과. 변경 Python 파일 py_compile 통과.

발견·수정한 오류:

1. Og 기본 Treutler 격자에서 원소별 xi 표 범위를 벗어나 IndexError. 반경 표 없이 생성하는 공식 Gauss-Chebyshev/Becke 격자를 도입했다. 없는 원자 반경 데이터를 임의 생성하지 않았다.
2. 격자 전자 수 함수 호출을 모듈 함수에서 NumInt 인스턴스 메서드로 수정했다.
3. GUI 주기계 지원표 키를 실제 카탈로그의 `valence_electrons`로 정정했다.
4. 지원 조합 총합 기대값의 산술 오타654를554로 정정했다. 개별 데이터 비교는 처음부터 일치했다.

### 경계 원소와 수치 격자

실제 Au⁺ PBE/def2-SVP/Newton, grid_level2, conv_tol1e-8: **수렴**, E=-135.389673188431 Ha,
전체78/코어60/명시18전자, Tr(DS)=18. `outputs/all_elements_validation/Au/results.json`.

최초 Og PBE/Dyall-v2z/scalar X2C/Gaussian/Newton:

| 방사/각도 격자 | 총에너지 [Ha] | 설명 |
| --- | ---: | --- |
| 90 / 302 (초기 level2) | -54565.771360841805 | 기존 PySCF 점 수, 격자 전자 수 오차6.75e-5 |
| 105 / 434 (초기 level3) | -54565.716064599197 | 에너지 차이 +0.055296242608 Ha, 격자 불충분 |
| 270 / 302 (최종 level2) | -54565.739878110828 | 방사3배, 격자 전자 수 오차 -6.68e-13 |
| 315 / 434 (최종 level3) | -54565.739872018950 | 방사3배, 격자 전자 수 오차 -2.42e-13 |

전자 수가 맞고 SCF가 수렴해도 에너지 격자 수렴이 보장되지 않는 사례다. 이 결과에 따라 universal의
Z>36 방사 점 수를 각 수준 기본값의 3배로 늘렸다. 낮은 격자 결과를 정확한 연구 기준값으로 사용하지 않는다.
270점 실행: `outputs/all_elements_validation/og_dense_grid/20261028_233425_o1rpk6lk`.
315점 실행: `outputs/all_elements_validation/og_dense_grid/20261028_233536_a512nxp9`.
최종 두 수준의 차이는 **6.09187555e-6 Ha (약0.000166 eV)**다. 이 비교는 해당 Og 계산의 수치 점검이며
기저·SOC·상관 근사 또는 다른 원소/분자의 정확도 검증이 아니다. 두 실행의 소요 시간68.27/114.55초는
다른 검증 작업과 동시 실행한 값이므로 성능 벤치마크로 사용하지 않는다.

최종 입력을 재현하는 스크립트는 `scripts/validate_heavy_elements.py`다. 각 실행에 고유 폴더와
입력·로그·전체 전자구조·결과 보고서를 남긴다. SOC, 이전자 picture-change, 바닥상태 자동 탐색은 미구현이며
실험/고정밀 상대론 기준값 대비 정확도는 이번 검증 범위에 포함하지 않는다.

최종 코드의 경계 원소 실행은 모두 수렴했고 전자 수·스핀·MO 직교성·ECP 코어 검사도 통과했다.
보고서: `outputs/all_elements_validation/20261028_233531_364618/report.json`.
공통 PBE/Newton, conv_tol1e-8, 4스레드, 메모리 목표4000MB, CUBE 축당16점. 스핀은 모두0이다.

| 계 | 기저·핵·상대론 | grid_level | 명시적 전자 | 총에너지 [Ha] |
| --- | --- | ---: | ---: | ---: |
| Cs⁺ | def2-SVP / 코어46 ECP | 2 | 8 | -20.000996063261 |
| Ce⁴⁺ | Dyall-v2z / Gaussian / scalar X2C | 2 | 54 | -8856.572163926809 |
| U⁶⁺ | Dyall-v2z / Gaussian / scalar X2C | 2 | 86 | -27982.333638639055 |
| Og | Dyall-v2z / Gaussian / scalar X2C | 3 | 118 | -54565.739872019010 |

Ce/U/Og는 보강된 universal 격자를 사용했다. 전자 수 Tr(DS)는 정수와 약1e-14 이내에서 일치했다.
이온은 f-block 전체전자 경로와 원소별 ECP를 검사하기 위해 선택했으며, 중성 원자의 모든 스핀 상태를
비교했다는 뜻은 아니다. 이전 낮은 방사 격자의 보고서도 삭제하지 않고 조건별로 보존했다.

## 0.6: 주기계와 MLIP 연결

- `test_periodic`: 입력/선정/구조 지문/미수렴 제외 6개, 실제 PySCF 엔진 2개 통과.
  Cu/Ag/Au/Pd/Pt의 기저·GTH 원자가 전자 수와 k점 계산의 **셀당 스핀** 보존을 검사했습니다.
- `test_periodic_gui`: 고배율 화면 입력 검사와 Windows→WSL 주기계 계산·그래프 표시 2개 통과.
- 기존 분자 입력 7개를 다시 실행해 통과했습니다. 아래 0.5의 전체 물리 검증 기록은 별도로 보존합니다.
- 실제 Si diamond 2원자, 2×2×2 k점: SCF 수렴, Tr(DS)=8.
  `outputs/periodic_validation/silicon_bulk/results.json`.
- 실제 Pt(111) 2층 2원자, 2×2×1 k점, 2D Coulomb: SCF 수렴, Tr(DS)=36.
  `outputs/periodic/20261028_135222_ztykquao/results.json`.
  최초 계산은 최종 수렴 검사 실패로 미수렴을 기록했고, 같은 허용오차의 제한된 재시도 경로를 검증했습니다.
- MACE-MP medium, float64, RTX 5080 Laptop: Ag8Au7Cu7Pd7Pt7+H 37원자, 36개 시작 후보 모두 H 완화 수렴.
  기본 e3nn 탐색 206.64초, 모델 로딩 34.88초 별도. 고정 표면입니다.
  `outputs/surfaces/hea_seed0_validation/screen_20261028_182228_578845/screen.json`.
  0.15 Å 기준으로 최종 H 위치를 묶으면 서로 다른 위치 19개입니다.
  가속을 켠 전체 탐색도 36개 모두 수렴, 같은 19개 위치로 분류됐고 선택 후보도 일치했습니다.
  탐색 97.78초, 모델 로딩 33.08초 별도.
  `outputs/surfaces/hea_seed0_validation/screen_20261028_183335_501532/screen.json`.
  전체 탐색 두 실행의 동시 CPU 부하 등은 통제하지 않았으므로, 가속비는 아래 별도 반복 측정을 기준으로 합니다.
- 작은 Pt2+H의 MLIP 상위 fcc/hcp 후보 두 개를 Windows→WSL로 실제 PBE 단일점 검증했고 둘 다 수렴했습니다.
  `outputs/surfaces/pt_small_validation/screen_20261028_181948_226356/dft_20261028_182225_448421/comparison.json`.
  MLIP와 DFT의 순위가 달랐습니다. 이 계산은 Gamma / 적분 수준 1의 **기능 검사**이므로,
  물리적 참값이나 MLIP 정확도 성적표로 사용하면 안 됩니다. 기저·k점·두께 수렴은 미검증입니다.
- cuEquivariance: 같은 모델 파일 SHA256·float64·GPU와 37원자 구조 3개, 각 7회 반복, ASE 캐시 초기화,
  CUDA 동기화, 각 구조 워밍업 후 에너지+힘을 측정했습니다. 중앙값 87.58 ms → 33.67 ms (**2.60배**).
  최대 에너지 차이 2.8422e-14 eV, 힘 성분 차이 1.3323e-15 eV/Å.
  `outputs/mace_acceleration_benchmark.json` 및 `benchmark_mace.py`로 재현합니다.

아직 검증하지 않은 것: 37원자 HEA 전체의 DFT 기준값, 물리적 흡착에너지 정확도,
SevenNet/CHGNet/MatterSim 실행·속도·정확도 비교, DFT 힘/응력/구조 최적화.
주기계 지원 추가와 같은 MACE의 GPU 가속을 확인했으며, 네 MLIP 전체를 능가했다는 결과는 없습니다.

## 0.5: 분자 DFT 검증

Windows Python 3.13 GUI와 Ubuntu 24.04 WSL Python 3.12.3 엔진을 연결해 검증했습니다.
엔진: PySCF 2.14.0, NumPy 2.5.4, Libxc 7.0.0.

## 자동 검증

| 묶음 | 테스트 | 결과 |
| --- | --- | --- |
| 기존 수소 물리/저장 | test_solver, test_workflow · 12개 | 이전 검증 통과 (0.4, 이번 변경 대상 아님) |
| 기존 수소 GUI | test_gui · 2개 | 이전 검증 통과 (2026-10-27, 이번 변경 대상 아님) |
| 분자 입력/물리 | test_molecular · 18개, WSL PySCF에서 실행 | 통과 |
| 분자 GUI | test_molecular_gui · 6개, 실제 Windows→WSL 실행 포함 | 통과 |

분자 물리 검증: H₂ 기준 에너지, He/H/H₂⁺의 전자 수와 스핀, MO의 S-직교성,
MO 계수/점유수로 밀도행렬 재구성, 기저함수·XC 교체, 물 분자의 Å↔bohr/평행이동 불변성.
SCF 최대 반복을 1로 제한한 물 분자(DIIS)와 Co²⁺(Newton)는 모두 `not_converged`로 저장됨을 확인했습니다.
GUI 검증은 좌표/XYZ 입력, JSON 저장, 수렴·미수렴 표시, 그림/ZIP 내보내기,
기록 복원, 엔진 오류 후 버튼 복구, 코발트 Newton 계산과 지르코늄 ECP 계산을 포함합니다.

지원 범위는 H~Zr (Z=1~40)이며 4주기 전체와 Rb·Sr·Y·Zr까지 포함합니다.
기저함수 7개 × 40개 원소의 280개 조합 중 241개는 내장 기저 데이터가 있고,
나머지 39개 조합은 입력 단계에서 차단됨을 확인했습니다.
Rb~Zr × def2-SVP/TZVP의 8개 조합은 대응 ECP의 코어 전자 수가 원자당 28개인지 확인했습니다.
혼합 분자의 기저 선택도 구성 원소의 공통 지원 범위로 제한합니다. 자동 대체하지 않습니다.
교환·상관 함수 4개의 Libxc 파싱은 이전 검증에서 확인했습니다.
이 확인은 각 조합의 과학적 정확도나 모든 금속 착물의 수렴을 보증하지 않습니다.

기존 확장 계산 검증: K·Ca 원자, Fe³⁺(2S=5)·Co²⁺(2S=3) 이온, KCl 분자.
각 결과의 MO 직교성·밀도행렬 재구성·전자 수·Nα−Nβ와 파일 저장을 확인했습니다.
Co²⁺는 기본 DIIS에서 미수렴했으며, 선택 가능한 2차 Newton 방식으로 수렴했습니다.
Newton의 중복 최종 콜백을 한 번만 기록하고, 제공되지 않는 밀도 변화 잔차는 null로 저장합니다.
미수렴은 반복 방식을 변경해도 성공으로 처리하지 않습니다.

0.5에서 추가한 13종은 Ni²⁺(2S=2), Cu⁺, Zn²⁺, Ga³⁺, Ge⁴⁺, As³⁺, Se⁴⁺, Br⁻, Kr,
Rb(2S=1), Sr, Y³⁺, Zr⁴⁺로 계산했고 모두 수렴했습니다. 괄호가 없는 상태는 2S=0입니다.
이는 선택한 상태의 수치 검증이며 모든 중성 원자·분자에 대한 바닥상태 정확도 검증은 아닙니다.

RbCl은 ECP 원소와 전체전자 원소가 섞인 분자로 검증했습니다. 전체 54전자 중 Rb 코어 28개를
제외한 명시적 전자 26개와 스핀 밀도를 검사했습니다. Newton 최대 1회로 제한한 RbCl도
`not_converged`로 기록됨을 확인했습니다.
ECP 계산의 CUBE 주석, 실제 원자번호, NPZ의 실제 원자번호/유효 핵전하/코어 전자 수,
그림의 Zr 원소 라벨 및 ECP 표시를 확인했습니다. ECP none과 Rb~Zr def2 기저의 잘못된 조합은 차단합니다.

## 실제 계산 예

| 계산 | 총에너지 [Ha] | 상태 |
| --- | --- | --- |
| H₂, 0.74 Å, PBE/STO-3G | −1.152072795159 | 수렴 |
| examples/water.xyz, PBE/def2-SVP | −76.272014472667 | 수렴, SCF 7회 |
| examples/co2plus.xyz, 전하 +2, 2S=3, PBE/def2-SVP, Newton | −1367.168991744082 | 수렴, 외부 반복 4회 |
| examples/zr4plus.xyz, 전하 +4, 2S=0, PBE/def2-SVP + def2-ECP | −44.023870246013 | 수렴, DIIS 5회 |

물의 `Tr(D S)`는 10.000000000000012로 입력 전자 10개와 일치했습니다.
Co²⁺의 `Tr(D S)`는 24.999999999999986, 밀도행렬의 Nα−Nβ는 3.000000000000001로
입력값 25전자와 2S=3에 일치했습니다. UKS 오비탈의 ⟨S²⟩는 3.7499999999999964입니다.
코발트 기본 정밀도(격자 수준 3) 실행은 Windows CLI→WSL 경로로 확인했고,
`outputs/molecular/20261028_132034_mfim2a4b`에 전체 결과를 보관했습니다.
Zr⁴⁺의 전체 전자는 36개, ECP 코어는 28개, 명시적 전자는 8개이며
`Tr(D S)`는 7.999999999999994입니다. 원자번호 40과 유효 핵전하 12를 구분해 저장했습니다.
기본 정밀도(격자 수준 3)의 Windows CLI→WSL 실행 결과는
`outputs/molecular/20261028_133015_ao7glxa3`에 보관했습니다.
ECP 계산과 전체전자 계산의 절대 총에너지는 직접 비교할 수 없습니다.
위 에너지는 명시된 좌표·기저·XC의 계산 결과이며 실험값이나 정확한 해가 아닙니다.
전체 재현 설정, 원시 결과와 엔진 버전은 `outputs/molecular`의 각 실행 폴더에 보관됩니다.

실행 명령과 테스트 범위는 [README_KO.md](../README_KO.md)를 참고하세요.
