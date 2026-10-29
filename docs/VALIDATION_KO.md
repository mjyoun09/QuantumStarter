# DFT 검증 기록 · 2026-10-29

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
