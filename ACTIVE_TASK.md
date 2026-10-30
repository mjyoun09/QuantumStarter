# 현재 작업

## 작업 식별 정보

- 작업 ID: QS-20261029-NATIVE-DFT
- 작업명: 첨부 범용 DFT 요구사항의 누락 구현, 자체 SCF 기반부터 단계별 검증
- 담당 AI: Claude (Claude Code, Opus 5.5) — 2026-10-29 20:14 사용자 요청으로 ChatGPT Codex의 HANDOFF_READY 인계서를 인수. 인수 시점 Codex 파일 수정 없음(마지막 기록 20:11), 실행 중 계산 없음 확인.
- 상태: HANDOFF_READY — 단계 9(0.9.0) 부분 완료: 자체 Newton SCF·주기계 밴드/DOS 구현·검증, **Cr·Co SCF 미해결**, 주기계 힘·응력은 PySCF GDF 미지원으로 보류. 다음은 아래 단계 10 재개 지점. 실행 중 계산·미저장 코드 없음. **GUI 시험은 QUANTUMSTARTER_TK_HIDDEN=1로만.**

### 단계 2: 원소별 기저 선택 (Claude, 완료 — CHANGELOG 0.8.1)

- 목적: 한 계산 안에서 원소마다 다른 기저를 지정(예: O=def2-tzvp, H=def2-svp). 인계 "정확한 다음 구현 지점" 1번.
- 설계: `DFTConfig.basis_overrides: {원소: 기저}` 추가. `basis`는 지정 없는 원소의 기본값으로 유지 → 기존 JSON/결과/체크포인트 호환. 같은 값 지정은 정규화로 제거. 구조에 없는 원소 지정은 허용하되 계산에 쓰이지 않음을 경고.
- ECP/X2C: 원소별 실제 기저로 판정(def2 무거운 원소만 def2-ECP). auto 상대론·핵·격자는 하나라도 Dyall이면 X2C/Gaussian/universal. ECP+X2C 혼합은 기존처럼 거부.
- 변경 예정: molecular.py(설정·검증), engines/pyscf_engine.py(dict 기저 전달, 결과 `basis_by_element`, 경고), scf_control.signature(지정 없으면 기존 해시 유지), energy_analysis(원소별 기저·ECP가 항마다 같은지 검사), main.py(`--basis-element EL=BASIS`), molecular_app.py(원소별 기저 입력), 테스트, docs.
- 완료 기준: 입력 검증 단위 테스트, PySCF 직접 구성과 AO 수/에너지 일치, 자체 엔진과 PySCF 엔진의 혼합 기저 에너지 일치(기존 비교 허용오차), 반응 분석의 원소별 기저 불일치 거부, 기존 회귀 테스트 통과, 문서 갱신.
- 기준점: Git 저장소를 새로 만들고 인수 직전 상태를 커밋 `dd2871e`로 보존(코드 변경 없음). outputs/는 .gitignore 대상이라 Git에 없음.
- 진행(20:40경): 구현 완료 — molecular.py(`basis_overrides`, `basis_by_element`, `basis_label`, `unused_basis_overrides`, `uses_relativistic_basis`, `parse_basis_overrides`/`format_basis_overrides`/`describe_bases`; 원소별 기저 데이터·ECP·auto 상대론/핵/격자), pyscf_engine(dict 기저, 결과 `basis_by_element`, 혼합/미사용 경고), scf_control.signature(사용된 지정만 해시; 없으면 기존 해시 그대로), energy_analysis(원소별 기저+ECP를 항 사이에서 대조, 기존 `ecp_policy` 키는 이로 대체), main.py `--basis-element`, GUI 입력란.
- 신규 test_basis_overrides.py 7개 + test_energy_analysis 1개 추가: WSL `python -m unittest test_basis_overrides test_energy_analysis -v` 10개 통과(3.5초). 혼합 기저 물: PySCF 엔진 = PySCF 직접 dict 기저 구성(<1e-8 Ha), 자체 엔진 차 <2e-7 Ha, AO 16개, AuH(Au def2-ECP+H STO-3G) 코어 60/명시 20전자 수렴.
- 회귀: WSL `test_native test_molecular test_periodic` 37개 통과(49.6초). Windows `QUANTUMSTARTER_GUI_DFT_TEST=1` `test_molecular_gui test_workflow test_basis_overrides test_energy_analysis` 24 통과/2 skip(PySCF 필요, WSL에서 통과), 54.8초. 버전 0.8.1로 올린 뒤 `test_workflow test_gui test_solver` 14개 통과.
- CLI 실계산 보존: outputs/basis_override_validation/20261029_202537_953szltv(pyscf), 20261029_202543_ocuvwdck(native), 차이 3.1e-11 Ha.
- 문서 갱신 완료: README_KO, docs/ELEMENTS·NATIVE_DFT·VALIDATION(0.8.1 절)·REQUIREMENTS_STATUS, PROJECT_STATE(사용자 목표 = DFT 프로그램 자체 제작, Git 상태), CHANGELOG. CLAUDE.md에 `@AGENTS.md` import 추가(Codex 문장은 유지).

### 단계 3: 분산 보정 (Claude, 완료 — CHANGELOG 0.8.2)

- 사용자 결정(2026-10-29): "더 정확한 것". 해석: 같은 모델이면 직접 구현과 라이브러리는 결과가 같고 정확도는 모델이 정한다 → **DFT-D4(EEQ 전하 의존 C6, BJ 감쇠 + ATM 3체)**를 Grimme 그룹 공식 구현 `dftd4` Python 패키지로 계산. D3 자체 구현은 하지 않음.
- 계획: (1) WSL .venv-dft에 dftd4 설치·API/원소 범위/매개변수 조회 확인 → (2) `DFTConfig.dispersion` none|d4, XC별 감쇠 매개변수(pbe/pbe0/b3lyp; LDA 거부), 원소 범위 검증 → (3) engines/dispersion.py: PySCF mol의 bohr 좌표·원자번호·총전하로 에너지, 결과에 성분 `dispersion`·모델·매개변수·버전 기록, 총에너지 = SCF + 분산 → (4) 체크포인트 서명에서 dispersion 제외(밀도 무관), 반응 분석은 분산 모델 일치 검사 → (5) GUI/CLI → (6) 검증: 이원자 장거리 −E·R⁶ → s6·C6 점근, 단위(Å/bohr) 불변, 병진·회전 불변, 전하 전달, 패키지 없을 때 명확한 오류·기본 기능 유지, 기존 회귀 → (7) 문서/CHANGELOG/커밋.
- 선택 의존성: 패키지 없이도 dispersion=none 경로는 그대로 동작해야 함. 설치 목록은 별도 requirements 파일로.
- 환경 변경: WSL `.venv-dft`에 `pip install dftd4`(4.3.1, cffi/pycparser 동반) 및 `tomli`(2.5.0; dftd4.parameters가 TOML 파서 요구) 설치. requirements-dispersion.txt 신설.
- 실측 사실: dftd4 지원 원소 Z=1~103, 112~118(110개, Rf~Rg 거부). 감쇠 매개변수 pbe/pbe0/b3lyp 존재, lda/svwn KeyError. `data.sqrt_z_r4_over_r2`는 원자번호 인덱스 배열(0번 더미). 기본 실공간 cutoff: 2체 60 bohr에서 쌍을 버림(59 포함/61 제외 확인) → 고립 분자에는 최대 원자간 거리+1 bohr로 확장(절단 없음), 결과에 cutoff 기록. 2체 D4-BJ 식 수동 재현이 라이브러리와 기계 정밀도 일치.
- 구현(20:55경): engines/dispersion.py(d4_energy), molecular.py(`dispersion` none|d4, D4_FUNCTIONALS/D4_ELEMENTS 검증), pyscf_engine(SCF 전 계산·SCF 후 합산, 성분 `dispersion`, `scf_energy_hartree`, `dispersion` 상세, 경고, assumptions), scf_control.signature(dispersion 제외), energy_analysis(분산 모델·매개변수·라이브러리 버전 일치), GUI 콤보·결과 표시, CLI `--dispersion`.
- test_dispersion.py 7개 WSL 통과(9.2초): 라이브러리 지원 범위/매개변수 계약, 2체 식 재현(R=3~200 bohr, 1e-10), 회전·병진 불변, 전하 전달(EEQ 합=총전하), 두 엔진 SCF 통합 + dftd4 자체 PySCF 어댑터와 1e-12 일치, Å/bohr 입력 일치, 라이브러리 부재 시 명시적 오류·none 경로 정상.
- 회귀: WSL 전체 54개 통과(첫 실행 1 실패 = test_basis_overrides 과거 서명 재현이 새 필드 dispersion을 빼지 않은 시험 누락, 시험 수정 후 통과). Windows 실계산 GUI 포함 44개 중 37 통과·7 skip(PySCF/dftd4 필요 6개는 WSL 통과, 기존 주기계 실계산 GUI 1개는 별도 환경변수 필요). 버전 0.8.2 변경 후 빠른 시험 재통과.
- CLI 실계산: outputs/dispersion_validation/ (PySCF+D4, 자체+D4, D4 없음). 두 엔진 차 2e-12 Ha, D4 −1.953413873068e-4 Ha. 실행 당시 version 필드는 0.8.1(문자열 변경 직전).
- 문서: README(분산 절·예제), START_HERE(선택 설치), docs/NATIVE_DFT·VALIDATION(0.8.2)·REQUIREMENTS_STATUS, PROJECT_STATE, CHANGELOG 갱신. **단계 3 완료.**

### 단계 4: DFT 원자 힘 (Claude, 완료 — CHANGELOG 0.8.3)

- 계획: `calculation='gradient'` 추가. 자체 엔진은 신규 engines/native_gradient.py에서 해석적 기울기를 직접 조립(1전자 기저·핵연산자 미분, J/K는 int2e_ip1 셸 단위 블록으로 메모리 제한, LDA/GGA XC는 AO 2차 도함수, 에너지가중 밀도 W=½DFD(RKS)/ΣD_sF_sD_s(UKS), 핵반발). 격자 가중치 미분 없음(PySCF 기본과 동일)을 결과에 명시. PySCF 엔진은 `mf.nuc_grad_method()`. D4 기울기는 dftd4 grad=True. SCF 미수렴이면 힘을 내지 않음.
- 검증 계획: (a) **고정 격자** 중심 유한차분(격자를 기준 구조에 고정하면 해석식이 정확한 미분이어야 함) RKS/UKS × LDA/PBE/PBE0, (b) 같은 격자의 PySCF 기울기와 비교, (c) 움직이는 격자 유한차분과의 차이(격자 응답 크기)를 격자 수준별로 기록, (d) D4 기울기 유한차분, (e) PySCF 경로 ECP 실행, (f) 회귀·GUI.
- 결과 형식: `gradient` 블록(Ha/bohr 기울기, eV/Å 힘, 최대/RMS, 합력, 성분별(자체), 드라이버, grid_response=false).
- 구현(21:15경): engines/native_gradient.py, engines/dispersion.py(gradient=True), pyscf_engine `nuclear_gradient`(자체/PySCF, D4 합산, 미수렴이면 힘 없음+경고), molecular.py `CALCULATIONS`, main.py `--calculation`, GUI `계산 종류` 콤보와 최대 힘 표시.
- 첫 비교에서 자체 해석 기울기 = PySCF 기울기(같은 격자) 1e-15 수준으로 일치(물 LDA/PBE/PBE0, OH UKS PBE).
- 오류·조치: test_gradient 첫 실행 1 실패 — OH UKS(격자 1, 엄격 수렴)에서 **자체 SCF가 300회 진동하며 미수렴**(ΔD ~2e-5 정체; 거의 축퇴한 π 궤도 점유 교대로 추정, 기울기 문제 아님). 시험 대상을 비축퇴 라디칼 NH2로 교체해 통과. **미해결 기록: 자체 DIIS의 근축퇴 열린 껍질 수렴 문제**(Codex의 OH 관찰과 같은 계열 가능성). 차후 레벨 시프트/점유 고정 등 검토 필요.
- test_gradient 6개 WSL 통과(11초). scripts/validate_gradient.py → outputs/gradient_validation/20261029_211817/report.json: 고정 격자 FD 최대차 1.2e-9~3.3e-9 Ha/bohr(5조건 × 전 좌표), 움직이는 격자와 차(격자 응답 크기) 수준1~4: 6.0e-4/4.4e-5/1.1e-5/1.1e-6, 자체 vs PySCF def2-SVP 8.8e-13/5.7e-13, PySCF Newton·sfX2C 기울기 지원 확인.
- 회귀: WSL 60개 통과(56.7초). Windows 50개 중 38 통과·12 skip(11개 WSL에서 통과, 기존 주기계 GUI 1개). 새 실계산 GUI 기울기 시험 통과(7.8초). 버전 0.8.3 변경 후 빠른 시험 재통과.
- CLI 실계산: outputs/gradient_validation/cli/ (물 PBE-D4/def2-SVP, 두 엔진 최대 힘 1.168150/1.168146 eV/Å, 기본 수렴 기준 차 ~1e-7 Ha/bohr). 실행 당시 version 필드 0.8.2.
- 문서 갱신: README(원자 힘 절), docs/NATIVE_DFT·VALIDATION(0.8.3)·REQUIREMENTS_STATUS, PROJECT_STATE(알려진 문제 포함), CHANGELOG. 커밋 c9719b1. **단계 4 완료.**
- (인수 당시의 단계 4 사전 계획 문단은 구현과 달라진 내용 — "calculation은 single_point만 허용" — 이 있어 삭제했다. 실제 구현은 위 기록과 커밋을 따른다.)

### 단계 5: SCF 수렴 수정 + 구조 최적화 (Claude, 완료 — CHANGELOG 0.8.4; 사용자 결정: 최적화기 geomeTRIC)

1. **[완료] 자체 SCF 수렴 정체 수정.** 진단(재현 조건: OH `O 0 0 0; H 0.1 0 .97`, spin 1, PBE, STO-3G, native, grid 1, conv_tol 1e-12, 밀도 1e-10, 잔차 1e-9): 점유 교대가 아니었다(연속 β 점유공간 겹침 0.9999999999, β HOMO-LUMO 간격 0.027 Ha 일정). **원인은 PulayDIIS의 계수 상한 30**: 수렴 근처에서 오차 벡터가 거의 평행해 계수가 커지는 것이 정상인데, 300회 중 289회 이력을 버려 사실상 DIIS 없이 느린 Roothaan 반복이 됨. 밀도 혼합 0·레벨 시프트 0.1/0.2·DIIS 16은 효과 없음. 수정: B 행렬 고유분해 의사역행렬(상대 1e-14 이하 고유값만 제거), 상한 제거, 비유한일 때만 리셋. 결과 OH 25회 수렴, H2O/NH2/O2 반복 수·에너지 불변. 회귀 시험 test_native `test_near_degenerate_open_shell_converges`(수정 전 코드에서 실패 확인).
   - **물리적 발견**: 수렴한 자체 OH 해는 PySCF 해보다 8.27e-6 Ha 높고, PySCF `stability()`로 **내부 불안정(안장점)** 판정. 불안정 방향을 따라가면 더 낮아짐. PySCF 자체 해도 불안정 판정 → 이 계는 다중 UKS 해. Codex 기록의 OH 2.4e-6 Ha 독립 차이도 같은 계열로 보인다. **수렴 ≠ 최소**; SCF 안정성 분석 기능은 미구현(향후).
   - 남은 한계: Fe 원자(STO-3G, 2S=4, PBE)는 수정 전후 모두 300회 미수렴(PySCF는 수렴).
2. 구조 최적화: **사용자 결정: geomeTRIC**(2026-10-29 "ㄱㄱ"). WSL .venv-dft에 geometric 1.1.1(+networkx 3.7) 설치, requirements-optimize.txt.
   - 구현(진행 중): engines/optimizer.py(geomeTRIC Engine.calc_new → calculate_molecule(gradient) 단계별 steps/step_NNN, 이전 단계 밀도를 초기 추정 `guess_density`로 전달, 단계 SCF 미수렴 시 중단, 최종 구조 재계산, optimization.xyz/optimized.xyz/final_settings.json, results.json `optimization` 블록, settings.json = 요청), engines/geometric_log.ini(로그는 파일만), pyscf_engine `guess_density`(initial_guess_details.resolved='previous_geometry_density', restarted는 체크포인트만), molecular.py `optimize`·`OPT_CONVERGENCE`(default=geomeTRIC 기본, tight=10배)·`opt_max_steps`, signature에서 opt 필드 제외, CLI `--calculation optimize --opt-convergence --opt-max-steps`, GUI 콤보/입력·지표.
   - 첫 실계산(tight, PBE/def2-SVP, grid 3): H2 자체 6단계 r=0.768334 Å, H2O 자체/PySCF 각 5단계 r(OH)=0.974677 Å ∠=102.044013°, **PySCF 자체 geomeTRIC 드라이버(pyscf.geomopt) 결과와 소수 6자리 일치**. 최종 재계산 에너지 차 ≤2.4e-10 Ha.
   - test_optimize.py 4개 통과(PySCF 자체 geomeTRIC 드라이버와 구조 비교, 기록, 단계 한도·SCF 실패·중단·라이브러리 부재, 작업자 경유). 실계산 GUI 최적화 시험 추가. 회귀: WSL 65 통과(첫 실행 1 실패 = test_gradient가 'optimize'를 잘못된 값 예시로 쓴 시험 문제 → 'hessian'으로 교체), Windows 56 중 41 통과·15 skip(14개 WSL 통과 + 기존 주기계 GUI 1).
   - CLI 실계산: outputs/optimization_validation/cli/ (왜곡 물 examples/water_distorted.xyz, PBE-D4/def2-SVP, 두 엔진 4단계, r(OH) 0.974723 Å, ∠ 102.0410°, 동일). 실행 당시 version 0.8.3.
   - 문서: README(구조 최적화 절), START_HERE(선택 설치 8), docs/NATIVE_DFT·VALIDATION(0.8.4)·REQUIREMENTS_STATUS, PROJECT_STATE, CHANGELOG. **단계 5 완료.**
   - 미구현으로 남긴 것: 최적화 체크포인트 재시작, 제약/스캔/TS, 진동수 확인.

### 단계 6: Hessian·진동수·ZPE (Claude, 완료 — CHANGELOG 0.8.5)

- 계획: `calculation=frequencies`. 해석적 힘의 중심 유한차분 Hessian(3N 좌표 × ±h, 기본 h=0.005 bohr, 각 변위 계산을 displacements/에 기록, 기준 밀도로 초기 추정, 미수렴 시 중단), 대칭화·비대칭 크기 기록, 최빈 동위원소 질량 가중, 병진·회전(선형 5개) 공간을 정확히 제거한 내부공간에서 대각화 → 3N−6(5) 진동수(cm⁻¹, 허수는 음수), 정규모드, ZPE, 진동 열역학(조화 근사, T 명시). 검증: PySCF 해석적 Hessian과 진동수 비교, 우리 Hessian을 PySCF `thermo.harmonic_analysis`에 넣어 투영·질량가중 독립 확인, H2 에너지 2차 미분으로 독립 확인, 평면 NH3 허수 모드, 선형 분자 모드 수.
- 발견: PySCF `atom_mass_list(isotope_avg=False)`는 질량수 정수(16, 1, 1)를 돌려준다 → 사용하지 않고 `pyscf.data.elements.COMMON_ISOTOPE_MASSES`(1.007825, 15.994915 등 정확한 동위원소 질량) 사용.
- 구현: engines/frequencies.py(`frequency_molecule`, `harmonic_analysis`, `vibrational_thermo`), molecular.py(`frequencies`, `freq_step_bohr`, `thermo_temperature_k`), signature 제외, dispatch, GUI 입력·지표, CLI `--calculation frequencies --freq-step --temperature`.
- 개발 비교(PBE/def2-SVP grid3, tight 최적화 구조): 자체 1608.59/3691.20/3790.51 vs PySCF 해석적 Hessian 1608.64/3690.72/3790.50 cm⁻¹(≤0.5); 우리 Hessian → PySCF harmonic_analysis 동일(소수 3자리); H2 4277.63 vs d²E/dr² 4277.72 cm⁻¹; 병진·회전 잔여 최대 ~18 cm⁻¹(격자 응답 미포함 탓, 정확히 제거됨); 물 계산 ~9초.
- 오류·조치: T=1 K에서 expm1 오버플로 → Cv가 inf/inf=NaN 가능(JSON 저장 실패 위험). q=e^(−x) 기반 안정식으로 수정, 경고를 오류로 둔 시험 통과.
- test_frequencies.py 6개 WSL 통과(25.8초).
- GUI 실계산 진동수 시험 추가·통과. 회귀 중 **기존 수소 GUI 시험 간헐 실패**(Windows 전체 실행 3회 중 2회, 110행: 기록 2번째 항목이 Z=2 실행) → 원인: 실행 폴더 이름 `YYYYmmdd_HHMMSS_무작위`를 이름 역순 정렬해 같은 1초 안의 두 실행 순서가 무작위(앱 기록 목록의 실제 버그). workflow.py·dft_workflow.py 폴더 접두사에 `%f`(마이크로초) 추가. 실패 조합 5회 반복 통과. 사용자 config.json은 변경되지 않았음을 확인.
- 최종 회귀: WSL 71 통과(82.9초), Windows 63 중 44 통과·19 skip(18개 WSL 통과 + 기존 주기계 GUI 1).
- CLI 실계산: outputs/frequency_validation/cli/(0.8.4 최적화 물 구조, PBE-D4/def2-SVP; 자체 1608.48/3690.73/3790.01, PySCF 1608.43/3690.50/3789.62 cm⁻¹, minimum, ZPE 54.366/54.362 kJ/mol). 실행 당시 version 0.8.4.
- 문서: README(진동수 절·예제), docs/VALIDATION(0.8.5)·REQUIREMENTS_STATUS, PROJECT_STATE, CHANGELOG. **단계 6 완료.**

### 단계 7: 전부 진행 (Claude, 완료 — B 0.8.6, A 0.8.7, C 0.8.8; NEB는 단계 8로)

- B-2/3) 자체 SCF 반복 수·전이금속 (커밋 예정 0.8.6 일부): 진단 — 기본 밀도 혼합 0.15가 쉬운 계에 3~4회, core 초기 추정이 2~3회 추가. 변경: (a) native `auto` 초기 추정 = PySCF `init_guess_by_minao`(원자 밀도 중첩, 외부 보조로 명시; 투영 밀도는 전자 수가 정확하지 않아 첫 Fock 생성에만 쓰고 그 점유 궤도에서 시작 — 처음엔 밀도를 직접 썼다가 취소 체크포인트의 전자 수 검사 실패를 test_native가 잡아 수정), (b) `native_damping_residual`(기본 1e-3): 교환자 노름이 이보다 클 때만 혼합. 결과 물 14→9, CH4 14→7, NH2 14→10, O2 10→8회; Fe 원자(STO-3G, 2S=4) 미수렴→수렴. 
  - 미해결(진단 기록): Cr(2S=6)·Co 원자는 어떤 혼합·추정에서도 정수 점유 미수렴. 페르미 스미어링(σ=0.01~0.03 Ha) 상태는 13~25회 수렴하지만 정수 점유로 되돌리면 다시 진동. 진동 에너지가 PySCF 수렴해보다 ~1e-5 Ha 낮음(PySCF 해도 최저가 아닐 수 있음). 축퇴 d 껍질의 aufbau 고정점 부재로 추정 → 2차(Newton) 해법 또는 분수 점유 모델 필요(미구현).
- B-1) SCF 안정성 분석 완료: engines/native_stability.py(실수 내부, check/follow), test_stability 3개 통과, PySCF 판정 일치, OH follow로 9.3e-7 Ha 낮은 안정해. WSL 75 통과, Windows 숨김 27 통과(2 skip). 문서·CHANGELOG 0.8.6. **7-B 완료.** 다음: 7-A 자유에너지(병진·회전·전자 열역학, PySCF hessian.thermo.thermo와 비교, 반응 분석 ZPE/G 옵션) → 7-C.
- A) 자유에너지 완료(0.8.7): engines/frequencies.py `ideal_gas_thermo`·`symmetry_number`, energy_analysis THERMO_FIELDS, 설정·GUI·CLI. PySCF thermo가 평균 질량을 쓰는 점 확인(비교 시 동일 질량 사용). test_frequencies 8 통과(경고=오류), WSL 77, Windows 숨김 31(+6 skip). **7-A 완료.** 다음: 7-C 제약·스캔(geomeTRIC constraints), NEB/TS.
- 0) 숨김 GUI 시험 모드 추가: test_molecular_gui `QUANTUMSTARTER_TK_HIDDEN=1`(루트 withdraw + Toplevel 생성 즉시 withdraw, 화면 배치 단언 2곳만 생략). 실계산 포함 14개 통과(79초). 수소·주기계 GUI 시험은 숨김 모드가 없어 실행하지 않음.

- (계획 당시 후보 A/B/C 목록은 아래 완료 기록으로 대체됨.)
- C) 제약·스캔·전이상태 완료(0.8.8): engines/optimizer.py(`opt_constraints`→constraints.txt, transition=True·hessian='first', `scan_molecule`), molecular.py `parse_scan`. 개발 오류: 사용자 제약이 내부 gradient 단계 설정으로 넘어가 검증에 걸림 → 내부·최종 설정에서 제약 제거. 첫 python3 실행이 stdin 대기로 멈춤(내 명령 실수, 파일 무변경 확인 후 .venv python 재실행). test_constraints 4 통과, WSL 81, Windows 숨김 27(+3 skip). **7-C 완료(NEB 제외).**

### 단계 8: NEB (Claude, 완료 — CHANGELOG 0.8.9)

- 구현: engines/neb.py(geomeTRIC ElasticBand/OptimizeChain/NEBParams, 선형 보간, 이미지별 밀도 재사용, neb/calc_NNNN 기록, 사이클 기록 save_callback, neb.csv·neb_path.xyz, 최고 이미지 재계산), 설정 neb_product·neb_images·neb_maxg_ev_per_angstrom, GUI·CLI.
- 확인한 사실: geomeTRIC NEB는 보간하지 않고 입력 프레임을 이미지로 쓴다; 양 끝 고정; maxg/avgg eV/Å; OptimizeChain은 maxcyc 도달 시에도 반환 → 반환 후 기준 직접 판정.
- 검증: test_neb 3 통과. NH₃ 반전 7 이미지 5회 수렴, 대칭 프로필(54.986 kJ/mol), 최고 이미지 평면, TS 정밀화 장벽 54.985 kJ/mol. WSL 84, Windows 숨김 32(+4 skip).
- 사고·조치: Windows Python write_text가 7개 소스를 CRLF로 바꿈 → 커밋 전 LF 복원(바이트 검사). 첫 CRLF 검사 명령은 오탐. PROJECT_STATE에 주의사항 기록.

### 단계 9: 자체 2차 SCF + 주기계 밴드·DOS (Claude, 부분 완료 — CHANGELOG 0.9.0)

- 9-1 계획: native 궤도 회전 신뢰영역 Newton(Steihaug CG, Hessian-벡터 곱은 native_stability의 기울기 중심차분 재사용). 점유는 aufbau를 강제하지 않고 현재 점유 궤도를 유지; 마지막 정준화는 점유/가상 블록 각각 대각화. aufbau 위반이면 경고. scf_solver=newton을 native에서 허용. 검증: Fe·Cr·Co 수렴, PySCF newton과 에너지 비교, 수렴 해의 안정성(최저 고유값 ≥ 0), 물은 DIIS와 같은 에너지.
- 9-1 진행(2026-10-30): engines/native_newton.py(Steihaug 신뢰영역 Newton, 점유는 궤도 번호로 유지, 수렴 시 최저 궤도 Hessian 고유값 검사 → 음수면 그 방향 선 탐색으로 안장점 탈출 후 계속, 블록 정준화, aufbau 위반 경고), native_scf 분기, molecular.py에서 native newton 허용, pyscf_engine의 PySCF newton 래퍼를 PySCF 엔진에만 적용(처음엔 native에도 적용되던 버그), 결과 `native_newton`. test_native의 'native+newton=오류' 단언 갱신. test_newton 2 통과(245초, 느림): 물 Newton=DIIS(1e-8), OH 안장점에서 시작해 더 낮은 안정 최소점.
  - 프로브: 물 9회(36.8초)·Fe 57회(183초) DIIS와 같은 에너지. Cr: 탈출 없을 때 안장점(최저 고유값 −0.42, E 0.13 Ha 높음)에 '수렴' → 탈출 추가 후 100회 미수렴이나 E −1032.808862(PySCF DIIS 수렴해 −1032.808849보다 낮음). 300회 실행은 백그라운드로 확인 중(scratchpad/cr300.log). PySCF newton 비교 프로브는 내 래퍼 사용 오류로 무효(물 −75.795) → PySCF DIIS 값을 기준으로 씀.
  - Cr 300회(scratchpad/cr300.log): 8회째 안장점 탈출(λ=−0.42) 후 회당 ~1e-12 Ha 감소, 기울기 4~7e-5 진동으로 미수렴(667초). 도달 E −1032.8088625(PySCF 수렴해보다 1.4e-5 낮음), aufbau 위반 없음. **Cr·Co 미해결 유지.**
- 9-2 진행: PySCF 2.14 실측 — GDF 주기계 기울기 없음(`'GDF' object has no attribute 'get_j_e1'`), 응력은 FFTDF 전용(AssertionError) → 힘·응력 미구현 결정(엔진은 slab 2D 쿨롱 때문에 GDF). get_bands는 GDF에서 동작. 구현: engines/periodic_analysis.py(DOS 항상, 밴드는 band_path 지정 시; ASE bandpath, slab은 2D), PeriodicConfig dos_sigma_ev·band_path·band_points, dos.csv·bands.csv, CLI --band-path/--band-points. Si 2×2×2 tight: DOS 적분 51.999999/52, 간접 갭 0.488 eV(Γ→Γ-X), 직접 2.51 eV, 메쉬 일관성 1.1 meV(conv_tol 1e-10에서도 → get_bands 경로 고유 차이로 판단). test_periodic_analysis 3 + test_periodic 9 통과. slab 밴드 경로는 실계산 미검증.
- 9-2 계획(원래): 주기계 밴드(고대칭 k 경로)·DOS(k 메쉬 고유값 Gaussian 확장) 후처리. PySCF PBC 기울기/응력 지원 여부를 실제로 확인 후 결정.

- 최종 회귀: WSL 89 통과(68+21), Windows 숨김 28 통과(+3 skip WSL 통과). 문서·CHANGELOG·PROJECT_STATE 0.9.0.

### 단계 10 재개 지점 (미착수)

1. Cr·Co 원자 자체 SCF: Newton이 평평한 방향에서 기울기 4~7e-5로 진동(신뢰영역 반경 진동 추정). 후보: 해석적 궤도 Hessian(XC 커널 f_xc 필요) 또는 FD 단계 크기·CG 정밀도 조정, 분수 점유(앙상블) 모델과의 비교. 재현: `DFTConfig(atoms=Cr 0 0 0, spin=6, xc=pbe, basis=sto-3g, engine=native, grid_level=3, scf_solver=newton, max_cycle=300)`.
2. Newton 속도: Hessian-벡터 곱 2회 Fock 생성 → 해석적 선형응답 또는 근사 대각 전처리.
3. 주기계 힘·응력: FFTDF 기반 별도 경로(bulk 한정) 검토 — 현재 GDF는 PySCF 2.14에서 미지원.
4. slab 밴드 2D 경로 실계산 검증, IRC, quasi-RRHO·용매, 연구용 HEA DFT 기준값.

(원래 단계 9 재개 목록)

1. IRC(TS→양쪽 최소점 연결 확인).
2. 자체 SCF: Cr·Co 원자 정수 점유 미수렴(2차/Newton 해법), 삼중항(RKS→UKS)·복소 불안정 검사.
3. 자유에너지: quasi-RRHO, 용매·표준상태(1 M) 보정.
4. 주기계: 힘/응력/셀 최적화, 밴드·DOS, D4 주기계(현재 분자만).
5. 연구용 HEA DFT 기준값과 MLIP 비교(PROJECT_STATE 향후 방향).
- 아래 "단계 8 재개 지점" 이후의 옛 문단은 계획 기록이다.

1. **NEB**(사용자 "전부 해"의 남은 항목): geomeTRIC 1.1.1 `geometric.neb` — `ElasticBand(molecule, engine, tmpdir, params, coords=None, plain=0)`, `OptimizeChain(chain, engine, params, save_callback=None)`, `NEBParams`. 우리 `optimize_molecule`의 내부 `Engine.calc_new`를 이미지 계산에 재사용 가능. 필요: 생성물 구조 입력(예: `neb_product_xyz` 또는 다중 프레임 XYZ), 이미지 수, 이미지별 기록(neb/image_NN/step_NNN), 최고 에너지 이미지→`transition_state` 정밀화 연결. 검증 후보: H 이동 반응 또는 NH₃ 반전(두 피라미드 구조 사이 → 최고점이 평면 TS 근처, 위 TS 에너지와 비교).
2. IRC(선택), 다차원 스캔(선택).
3. 미해결 기록: Cr·Co 원자 자체 SCF 정수 점유 미수렴(2차 해법 필요), 삼중항·복소 SCF 불안정 미검사, quasi-RRHO 없음.
4. 이후: 주기계 힘/응력/밴드/DOS → 연구용 HEA DFT 기준값.
- **검증 규칙(사용자 요청)**: GUI 시험은 `QUANTUMSTARTER_TK_HIDDEN=1`로만. test_gui/test_periodic_gui는 화면에 창을 띄우므로 실행하지 않는다.

### 단계 3 결정 전 기록 (참고용, 결정은 단계 3 절 참조)
- 결정 필요: (a) Grimme D3(BJ)를 자체 구현 — 매개변수 표(C6 참조값, r0, r2r4)를 공개 데이터에서 가져와 에너지/좌표수 식을 직접 작성, 사용자의 "DFT를 직접 만든다" 목표에 부합, 독립 기준값은 simple-dftd3로 검증. (b) `simple-dftd3`/`dftd4` Python 패키지를 .venv-dft에 설치해 호출(빠르지만 외부 경계 추가). 어느 쪽이든 결과에 별도 에너지 성분 `dispersion`과 모델·매개변수 출처를 기록하고, 패키지 없이도 기본 기능은 동작해야 한다.
- 영향 파일 예상: molecular.py(DFTConfig.dispersion 필드·검증), engines/pyscf_engine.py(성분 합산·출처), engines/native_scf.py(성분 표시만; SCF에는 영향 없음 — D3는 밀도 비의존), scf_control.signature(분산은 밀도에 영향 없으므로 서명에서 제외 검토), energy_analysis.compatibility(분산 모델 일치 검사), GUI/CLI, 테스트.
- 이후 순서(인계 원문 유지): 힘(유한차분 검증) → 구조 최적화/제약/스캔 → Hessian/진동/ZPE → NEB → 주기계 응력/밴드/DOS.

## 작업 목적

- 사용자 요청: 첨부 Pasted text.txt의 미구현 기능을 실제 구현하고 검증.
- 기존 PROJECT_STATE/AGENTS/소스 확인 완료, ACTIVE_TASK 없음, Git 저장소 아님.
- 핵심 누락: 현재는 PySCF 외부 SCF 드라이버만 사용. 자체 밀도 구성·Fock 조립·DIIS·SCF·에너지 합산 없음.
- 프롬프트 9절의 우선순위에 따라 자체 분자 DFT 기반을 먼저 완성한다. 외부 엔진을 자체 구현으로 설명하지 않는다.
- 첫 단계 완료 기준: 자체 RKS/UKS·LDA/GGA/지원 global hybrid, 명시적 에너지 성분과 수렴 기준, GUI/CLI 엔진 선택, 체크포인트/중단/재시작, He/H2/H2O 등 PySCF와 동일 조건 비교 및 전자 수·에너지 합·잔차·초기추정 검증.
- 후속 요청 항목(힘/구조 최적화/스캔/진동/반응·흡착 분석/NEB/주기계 응력·밴드 등)은 요구사항 추적 문서에 실제 구현 상태를 남기고 순차 확장. 전체 요구사항을 첫 단계 완료만으로 완료 처리하지 않는다.

## 작업 계획

1. 요구사항 추적 문서와 외부/자체 경계 확정.
2. native_scf 모듈: PySCF는 기저·AO/ERI 적분·적분 격자, Libxc는 XC 함수값/미분만 제공. 밀도·XC 행렬·J/K·고유문제·Pulay DIIS·에너지·수렴은 직접 구현.
3. 설정/공통 출력/GUI/CLI에 자체 엔진 추가. 기존 PySCF와118원소 지원 보존. 자체 첫 단계는 H~Kr 비상대론 전체전자·점핵, 자원 한도 명시.
4. 중단 신호, 원자적 체크포인트 저장, 호환성 검사 후 재시작, 결과 출처 기록.
5. 동일 조건 외부 비교, 수렴 불변량, 강제 미수렴·중단·재시작·GUI, 기존 회귀 검사.
6. 상태/변경 이력/검증 문서 갱신. 미완료 요구사항은 숨기지 않고 인계서 보존.

## 진행 상황

완료: 첨부/소스 조사, 요구 추적 문서, 자체 SCF·설정·출력·중단/체크포인트·GUI/CLI, 구조 측정·프리셋, 결과 비교·반응/흡착에너지 분석. 실제 비교5종과 GUI 중단/재시작 검증 결과 저장.
첫 단계 완료: 최종 회귀 검사·영구 문서 동기화까지 마침. 실행 중인 계산이나 미저장 코드 없음.
미수행: 원소별 기저, 분산 보정, DFT 힘/최적화/제약/스캔/진동/열역학/NEB, 주기계 응력/밴드/DOS 등 첨부5절 후속 확장. 첫 단계의 완료와 전체 요청의 완료를 구분하므로 이 문서를 삭제하지 않는다.

## 구현 및 실행 기록

- 조사: PROJECT_STATE.md, AGENTS.md, engines/pyscf_engine.py, engines/__init__.py, molecular.py, dft_workflow.py, dft_worker.py, main.py, molecular_app.py, surface_workflow.py.
- 현재 `.venv` Windows GUI / `.venv-dft` WSL PySCF2.14 / `.venv-mlip` GPU 환경 구분. 기존 데이터·노트북 삭제/수정 없음.
- 신규 engines/native_scf.py: 자체 밀도/gradient/XC 행렬/J/K/고유문제/DIIS/수렴/에너지 성분. scf_control.py: 원자적 NPZ 체크포인트, 물리 조건 지문, 전자 수/점유수 검사, 중단 신호.
- 변경 molecular.py, engines/pyscf_engine.py, dft_workflow.py, dft_worker.py, molecular_app.py, main.py. native/pyscf 분리, 외부 기존 설정은 pyscf 유지. 외부의 density-based grid pruning을 끄고 동일 고정 격자 비교(변경 영향 회귀 검사 예정).
- 실제 물 LDA/STO-3G native: 13회 수렴, -74.73210308754511 Ha. 성분 합 일치, outputs/native_validation/first_water.
- test_native 첫 실행: RKS 6조건 외부 비교와 RKS Fock 에너지 미분, 초기추정/격자/미수렴, 중단·재시작 검사 통과. UKS 7서브케이스 오류: NumPy einsum의 ellipsis 축 합산 표현 오류. 명시적 np.sum(dm * K.T)로 수정, 재검증 진행.
- Git 상태: not a git repository. 데이터 삭제 없음.
- UKS 합산 수정 후 `test_native test_molecular test_periodic` 35개 모두 통과(40.237초). 자체9조건 기준 비교와 RKS/UKS×4XC 에너지 유한미분도 통과. Windows GUI 입력4개 통과, 기존 선택 실계산4개는 해당 실행에서 skip.
- geometry_tools.py 및 GUI: 구조 미리보기, 거리/각도, 실제 수치 변경을 표시하는 빠른/정밀 프리셋 추가. 새 테스트 및 검증 스크립트 scripts/validate_native.py 추가.
- 추가검증: 외부 DIIS/Newton 체크포인트 재시작과 입력 테스트3개 통과. 영구 비교 He/H2/H2O 에너지 차이 각각4.44e-16/4.44e-16/4.26e-14 Ha. OH 2스레드 독립 SCF는 둘 다 수렴하나2.4055e-6 Ha 차이로 엄격 비교 실패. 같은 밀도에서 Fock/에너지 평가를 비교해 구현 오류인지 퇴화 궤도/유한각도격자 상태 차이인지 조사 중. 보고서 경로 outputs/native_validation/20261029_194600_424557.
- GUI 실제 중단 테스트: 작은 STO-3G 계산이 UI 큐 처리 중 이미 끝나 수렴 결과를 반환함. 매 메시지마다 그래프를 그리던 루프를30ms 예산·배치당1회 갱신으로 수정, 충분히 긴 실제 계산으로 사용자 중단 검사 재실행 예정. 중단 로직 단위 실계산 검사는 이미 통과.
- 동일 OH 밀도에서 에너지 차이4.26e-14 Ha, Fock 최대차3.55e-15 Ha. 외부를 자체 최종 밀도에서 시작하면 같은 에너지에 수렴. 독립 해 차이는 별도로 보존하며 정확도/바닥상태 보증으로 일반화하지 않음.
- 체크포인트 schema2: 동일 물리·기저·격자 설정이면 native/pyscf 간 밀도 전달 허용, 원본 source_engine 기록. 엔진 자체 선택 대신 실제 물리 설정과 밀도 점유수·전자 수를 검사. 오래된 개발중 schema1은 버전 오류로 거부.
- GUI 체크포인트 생성 직후 실제 중단 버튼→재시작 통합 검사1개 통과(22.361초), 구조/프리셋 GUI1개 통과. 새 영구 비교5종 모두 완료: outputs/native_validation/20261029_195245_027252/comparison.json. OH 독립차2.40551e-6 Ha와 동일밀도재시작차2.84e-14 Ha 구분.
- energy_analysis.py / test_energy_analysis.py: 조건·원소/전하 수지·수렴 확인 후 반응/흡착에너지 합산, 참조파일SHA256·출처 보존. GUI 기록 비교/흡착 분석과 CLI 추가. 테스트2개 통과. 실제 H2−2H 분석 및 최종 전체 회귀·문서화 남음.
- 실제 H2−2H 분석 완료: -6.799190445 eV, outputs/analysis/20261029_195626_hr778v3r/analysis.json. LDA/STO-3G 합산·출처 검사이며 실험 결합에너지 정확도 검증이 아님.
- Windows CLI→WSL 자체 물 예제 완료: -74.7321051890996 Ha, 13회, outputs/molecular/20261029_195734_a6kkyr2a. 리다이렉트 출력의 한글 인코딩 문제를 CLI UTF-8 설정으로 수정.
- 최종 canonical 밀도와 이전 SCF 밀도 사이에도 에너지·밀도 수렴 기준을 모두 적용하도록 강화. 결과 이력과 최종 에너지·밀도의 일치 검사를 추가.
- 세션 출력이 잘려 이전 최종검사의 결과를 확인하지 못했으므로 통과로 간주하지 않음. 최종 수정 상태에서 물리/분석/기존 분자·주기계39개와 GUI10개를 재실행하고 outputs/native_validation/final_physics_tests.log, final_gui_tests.log에 원본 출력을 보존 중.
- 위 최종 검사가 실제 완료됨: WSL39개 모두 통과(66.843초), 실계산 GUI10개 모두 통과(64.352초, skip 없음), Windows test_workflow5개 통과(0.104초). 후속 코드는 CLI의 오류 안내만 수정했으며11개 변경 Python 파일 컴파일과 두 CLI UTF-8 한글 help 출력 검사를 통과했다.
- UTF-8 확인용 임시 명령의 상대 경로 `.venv/Scripts/python.exe`를 Windows subprocess가 찾지 못한 WinError2 발생. 앱 오류가 아니라 검사 명령의 경로 문제이며 Path.resolve()로 절대 경로를 사용한 재검사가 통과했다. 처음 실패를 CLI 검증 통과로 기록하지 않음.
- README_KO.md, START_HERE_KO.md, docs/NATIVE_DFT_KO.md, docs/ENERGY_ANALYSIS_KO.md, docs/VALIDATION_KO.md, docs/REQUIREMENTS_STATUS_KO.md, PROJECT_STATE.md, CHANGELOG.md를 실제 기능/한계/검증과 동기화. CHANGELOG에는 첫 단계 완료만 기록하고 전체 요청 완료는 주장하지 않는다.

## 다음 AI를 위한 인계

이 담당 AI의 첫 단계 작업은 저장·검증을 마쳐 인계 가능하다. 다른 세션의 동시 수정 여부는 별도로 확인한다. 새 담당자는 PROJECT_STATE → 이 문서 전체 → 실제 코드 순으로 확인하고 담당자·상태를 갱신한다. 자체 SCF를 처음부터 다시 구현하지 않는다.

우선 확인 명령(프로젝트 루트 PowerShell):

```powershell
Get-Content docs/REQUIREMENTS_STATUS_KO.md -Encoding UTF8
Get-Content outputs/native_validation/final_physics_tests.log -Tail 8
Get-Content outputs/native_validation/final_gui_tests.log -Tail 8
# 변경 후 필요한 범위의 검사를 선택한다. 기존54개를 이유 없이 반복할 필요는 없다.
wsl --exec .venv-dft/bin/python -m unittest test_native test_energy_analysis test_molecular test_periodic -v
$env:QUANTUMSTARTER_GUI_DFT_TEST='1'
.\.venv\Scripts\python.exe -m unittest test_molecular_gui -v
```

정확한 다음 구현 지점:

1. 첨부5절 순서의 남은 분자 확장: 원소별 기저 선택을 먼저 설계한다. molecular.py의 DFTConfig/기저 지원 검증, engines/pyscf_engine.py의 분자 구성, molecular_app.py의 기저 선택, scf_control.signature 및 energy_analysis.compatibility를 함께 검토한다. 기존 단일 문자열 JSON을 깨지 않으며 실제 원소별 기저/ECP가 다른 에너지는 섞지 않는다.
2. 분산력은 선택 의존성과 물리 모델(D3/D4 등)을 명시해 별도 에너지 성분으로 추가한 뒤 독립 기준값과 검증한다. 아직 라이브러리 설치나 모델 선택을 하지 않았다. 기본 기능이 선택 패키지 없이 동작해야 한다.
3. 힘/최적화 단계는 자체와 외부의 미분 구현 출처를 구분하고 먼저 에너지 유한차분으로 힘을 검증한다. 구조 제약·최적화 종료 시 힘 기준과 미수렴 구분을 구현한 뒤 스캔/Hessian/진동/ZPE로 확장한다. 현재 DFTConfig.calculation은 single_point만 허용한다.
4. 이후 주기계 응력/셀 최적화/밴드/DOS와 연구용 HEA DFT 기준값은 별도 작업 단위로 진행한다. 37원자 HEA의 독립 DFT 기준과 네 MLIP 우위 검증은 아직 없다.

주의·미검증 가정:

- 자체 엔진에서 PySCF RKS/UKS.kernel/get_veff/DIIS를 호출하면 핵심 요구를 만족하지 못한다. 적분/AO/격자와 Libxc 함수가 외부 경계다.
- 자체 H~Kr/160 AO 제한과 외부118원소 지원을 구분한다. 다중 전자 DFT 오차와 수소 방사방정식 정확해를 혼동하지 않는다.
- OH 독립 SCF의2.40551e-6 Ha 차이 원인은 완전히 확정하지 못했다. 같은 밀도 에너지/Fock 일치 및 외부 재시작 일치는 확인했다. 바닥상태 보증으로 쓰지 않는다.
- 체크포인트 schema2는 동일 물리 조건에서 엔진 간 전달을 허용하지만 DIIS 이력은 복원하지 않는다. 반응에너지 분석은 엔진 혼합을 거부하며 주기계는 같은 셀·k점만 허용한다.
- Git 저장소가 없어 커밋하지 않았다. 기존 사용자 연구 데이터와 노트북을 삭제/수정하지 않았다. 후속 작업이 남았으므로 ACTIVE_TASK.md를 보존한다.
