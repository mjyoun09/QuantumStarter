# 변경 이력

완료된 작업만 기록한다. 진행 작업은 ACTIVE_TASK.md, 장기 상태는 PROJECT_STATE.md를 참조한다.
규칙 등록 이전의 검증 이력은 [docs/VALIDATION_KO.md](docs/VALIDATION_KO.md)에 보존되어 있다.

## QS-20261029-NATIVE-DFT · 단계 9 (자체 2차 SCF + 주기계 밴드·DOS) 2026-10-30 · 0.9.0 · Claude

- 자체 엔진 신뢰영역 Newton SCF(engines/native_newton.py): Steihaug CG, 점유를 궤도 번호로 유지, 수렴 시 최저 궤도 Hessian 고유값 검사 후 음수면 선 탐색으로 탈출, 블록 정준화, aufbau 위반 경고, 결과 `native_newton`. 자체 엔진에서 scf_solver=newton 허용(PySCF 래퍼는 PySCF 엔진에만).
- 주기계 DOS(항상)·밴드(band_path 지정 시) 후처리(engines/periodic_analysis.py), dos.csv·bands.csv, PeriodicConfig dos_sigma_ev·band_path·band_points, CLI --band-path/--band-points.
- 주요 파일: engines/native_newton.py, engines/native_scf.py, engines/pyscf_engine.py, engines/periodic.py, engines/periodic_analysis.py, molecular.py, periodic.py, main.py, workflow.py, 신규 test_newton.py·test_periodic_analysis.py, test_native.py, README_KO.md, docs/*.
- 검증: Newton=DIIS(물), OH 안장점 탈출, Fe 수렴; Si DOS 적분·밴드 일관성(1.1 meV). WSL 89개, Windows 숨김 28개(+3 skip WSL 통과).
- **미해결**: Cr·Co 원자 SCF(Newton 300회도 미수렴), 주기계 힘·응력(PySCF GDF 미지원), slab 밴드 실계산 미검증.

## QS-20261029-NATIVE-DFT · 단계 8 (NEB) 완료 2026-10-30 · 0.8.9 · Claude

- `calculation=neb`: 신규 engines/neb.py — geomeTRIC CI-NEB(ElasticBand/OptimizeChain/NEBParams)에 자체/PySCF 엔진의 에너지·해석적 힘을 공급. 선형 보간 초기 경로(정렬 없음), 양 끝 고정, 이미지별 밀도 재사용, 모든 기울기 계산 기록(neb/calc_NNNN), 사이클별 에너지·힘 기록, neb.csv·neb_path.xyz, 최고 이미지 전체 재계산. 설정 `neb_product`·`neb_images`·`neb_maxg_ev_per_angstrom`, GUI·CLI(`--neb-product`, `--neb-images`).
- 주요 파일: engines/neb.py, engines/__init__.py, molecular.py, scf_control.py, main.py, molecular_app.py, workflow.py, 신규 test_neb.py, test_basis_overrides.py, README_KO.md, docs/VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: NH₃ 반전 대칭 경로, 최고 이미지 평면, TS 탐색 장벽과 0.001 kJ/mol 차. WSL 84개, Windows 숨김 32개(+4 skip WSL 통과).
- 한계: 선형 보간 초기 경로만, IRC 없음.

## QS-20261029-NATIVE-DFT · 단계 7-C (제약·스캔·전이상태) 완료 2026-10-29 · 0.8.8 · Claude

- 제약 최적화 `opt_constraints`(geomeTRIC `$freeze`/`$set` 문법, 실행 전 구문 검사, constraints.txt 기록), 이완 스캔 `calculation=scan` + `scan_spec`(distance/angle/dihedral, 점별 제약 최적화, scan.csv·scan.xyz, 상단 폴더에 마지막 점 출력), 전이상태 탐색 `calculation=transition_state`(geomeTRIC P-RFO, 초기 Hessian은 우리 힘 유한차분, 단계 종류 기록). GUI·CLI(`--constraints`, `--scan`).
- 주요 파일: engines/optimizer.py, engines/__init__.py, molecular.py(`parse_scan`, 필드·검증), scf_control.py, main.py, molecular_app.py, workflow.py, 신규 test_constraints.py, test_basis_overrides.py, README_KO.md, docs/VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: 각도 제약 100° 정확, 원자 고정 3.4e-5 Å 이내, H₂ 스캔 = 단일점 에너지(2e-10 Ha), NH₃ 반전 TS 허수 진동 정확히 1개. WSL 81개, Windows 숨김 27개(+3 skip WSL 통과).
- 미착수: NEB·IRC(geomeTRIC API 확인만), 다차원 스캔. ACTIVE_TASK.md 유지.

## QS-20261029-NATIVE-DFT · 단계 7-A (이상기체 자유에너지) 완료 2026-10-29 · 0.8.7 · Claude

- 진동수 계산에 RRHO 이상기체 열역학 추가: 병진·회전(선형/비선형/원자)·전자(다중도)·진동 항, H/S/G, ZPE 보정 에너지. 결과 최상위 `zpe_corrected_energy_hartree`·`enthalpy_hartree`·`gibbs_free_energy_hartree`. 설정 `thermo_pressure_pa`, `symmetry_number`(0=PySCF 점군 판정), GUI·CLI(`--pressure`, `--symmetry-number`).
- 반응·흡착 분석에서 ZPE/H/G 선택 가능, 이때 온도·압력·모델 일치 검사. 한계 문구를 선택한 에너지에 맞게 표시.
- 주요 파일: engines/frequencies.py, energy_analysis.py, molecular.py, scf_control.py, main.py, molecular_app.py, workflow.py, test_frequencies.py, test_basis_overrides.py, README_KO.md, docs/VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: PySCF thermo와 항목별 일치(H·G 1e-10 상대), WSL 77개, Windows 숨김 31개(+6 skip WSL 통과).
- 한계: quasi-RRHO·비조화·용매·표준상태 변환 없음, 들뜬 전자 상태 무시.

## QS-20261029-NATIVE-DFT · 단계 7-B (자체 SCF 품질) 완료 2026-10-29 · 0.8.6 · Claude

- 자체 SCF `auto` 초기 추정을 원자 밀도 중첩(PySCF minao 보조, 첫 Fock 생성에만 사용)으로, 밀도 혼합을 잔차 기준(`native_damping_residual`, 기본 1e-3)으로 변경(커밋 d5f9399). 반복 수 물 14→9, CH4 14→7, NH2 14→10, O2 10→8, Fe 원자 수렴.
- 신규 engines/native_stability.py: 실수 내부 궤도 Hessian 안정성 분석(기울기 직접 구성, Hessian-벡터 곱 중심차분, Lanczos), `native_stability=check|follow`(GUI·CLI `--stability`), 결과 `scf_stability`, 불안정 경고. follow는 불안정 모드 선 탐색 후 SCF 재수행.
- 테스트 인프라: 사용자 요청으로 GUI 시험 숨김 모드 `QUANTUMSTARTER_TK_HIDDEN=1`(커밋 0d7dc2f). 이후 화면에 창을 띄우는 시험은 실행하지 않음.
- 주요 파일: engines/native_scf.py, engines/native_stability.py, engines/pyscf_engine.py, molecular.py, scf_control.py, main.py, molecular_app.py, workflow.py, test_native.py, test_stability.py, test_molecular_gui.py, test_basis_overrides.py, README_KO.md, PROJECT_STATE.md, docs/NATIVE_DFT_KO.md·VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: 안정성 기울기·Hessian이 에너지 차분과 일치, 물 안정·OH 불안정 판정이 PySCF와 일치, follow로 OH가 9.3e-7 Ha 낮은 안정해로 이동. WSL 75개, Windows 숨김 GUI 포함 27개 통과(2 skip WSL 통과).
- 미해결: Cr·Co 원자 정수 점유 미수렴(진단 기록), 삼중항·복소 불안정 미검사.

## QS-20261029-NATIVE-DFT · 단계 6 (조화 진동수·ZPE) 완료 2026-10-29 · 0.8.5 · Claude

- `calculation=frequencies`(GUI `계산 종류`·변위·온도 입력, CLI `--calculation frequencies --freq-step --temperature`): 해석적 힘의 중심 유한차분 Hessian(6N+1회 계산, displacements/에 기록, 비대칭 크기 기록), 최빈 동위원소 정확 질량, 병진·회전(선형 5) 정확 제거 후 3N−6(5) 진동수, 정규모드, 허수 진동 판정(|ν|>20 cm⁻¹), ZPE·조화 진동 열역학(U/S/Cv). 정지점이 아니면 경고.
- 주요 파일: 신규 engines/frequencies.py, test_frequencies.py; molecular.py, scf_control.py, engines/__init__.py, main.py, molecular_app.py, test_molecular_gui.py, test_basis_overrides.py, test_gradient.py, workflow.py·dft_workflow.py(실행 폴더 이름 마이크로초), README_KO.md, docs/VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: 투영·질량가중이 PySCF harmonic_analysis와 1e-4 cm⁻¹, PySCF 해석적 Hessian 진동수와 ≤0.5 cm⁻¹(def2-SVP)·2 cm⁻¹ 허용 시험, H₂ 에너지 곡률과 0.09 cm⁻¹, 평면 NH₃ 허수 우산 모드, 두 엔진 CLI 물 PBE-D4 0.4 cm⁻¹ 이내. WSL 71개, Windows 44개(+19 skip 사유 기록) 통과.
- 수정한 문제: 저온 열역학 오버플로(NaN 가능), 같은 1초 안의 실행이 기록 목록에서 무작위 순서로 정렬되던 기존 버그(수소 GUI 시험 간헐 실패의 원인).
- 한계: 조화 근사, 병진·회전·전자 열역학과 G 없음, 해석적 Hessian 없음(6N회 계산). ACTIVE_TASK.md 유지.

## QS-20261029-NATIVE-DFT · 단계 5 (DIIS 수정 + 구조 최적화) 완료 2026-10-29 · 0.8.4 · Claude

- 자체 DIIS: Pulay 계수 상한 30을 없애고 B 행렬 고유분해 의사역행렬로 풀도록 수정(커밋 de7d245). OH에서 300회 중 289회 이력을 버려 정체하던 문제 해결(→ 약 25회 수렴), 다른 계 불변. 수렴한 OH 해가 PySCF 안정성 분석상 불안정임을 기록(다중 UKS 해).
- `calculation=optimize`(GUI `계산 종류`, CLI `--calculation optimize --opt-convergence --opt-max-steps`): geomeTRIC(TRIC)이 단계를 제안하고 자체/PySCF 엔진이 에너지·해석적 힘(D4 포함)을 계산. 단계별 전체 기록 steps/step_NNN, 이전 구조 밀도를 초기 추정으로 사용, 단계 SCF 미수렴 시 중단, 최종 구조 재계산, optimization.xyz/optimized.xyz/final_settings.json/geometric.log, results.json `optimization` 블록. geomeTRIC은 선택 의존성(requirements-optimize.txt).
- 주요 파일: 신규 engines/optimizer.py, engines/geometric_log.ini, requirements-optimize.txt, test_optimize.py, examples/water_distorted.xyz; engines/native_scf.py, engines/pyscf_engine.py(guess_density), engines/__init__.py, molecular.py(OPT_CONVERGENCE·opt 필드), scf_control.py, main.py, molecular_app.py, workflow.py, test_native.py, test_gradient.py, test_basis_overrides.py, test_molecular_gui.py, README_KO.md, START_HERE_KO.md, docs/NATIVE_DFT_KO.md·VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: PySCF 자체 geomeTRIC 드라이버와 물·H2 최적 구조 소수 6자리 일치(테스트 허용 2e-4), CLI 물 PBE-D4 두 엔진 동일 구조, test_optimize 4개, WSL 65개, Windows 41개(+15 skip 사유 기록) 통과. 상세 docs/VALIDATION_KO.md 0.8.4.
- 한계: 국소 최소점 후보(진동수 확인 없음), 제약/스캔/전이상태·최적화 재시작 없음, SCF 안정성 분석 없음, Fe 원자 자체 엔진 미수렴. ACTIVE_TASK.md 유지.

## QS-20261029-NATIVE-DFT · 단계 4 (원자 힘) 완료 2026-10-29 · 0.8.3 · Claude

- 분자 계산에 `calculation=gradient`(GUI `계산 종류`, CLI `--calculation gradient`)를 추가했다. 입력 구조에서 총에너지의 핵 좌표 기울기와 힘을 계산하며 구조는 바꾸지 않는다.
- 자체 엔진: 신규 engines/native_gradient.py에서 해석적 RKS/UKS 기울기를 직접 조립(1전자 기저·핵 이동, J/K 2전자 미분 셸 단위, LDA/GGA XC, 에너지가중 밀도 Pulay 항, 핵반발). PySCF 엔진은 nuc_grad_method. D4를 켜면 dftd4 기울기를 합산. SCF 미수렴이면 힘을 내지 않는다. 격자 응답 항은 포함하지 않으며 결과·경고에 명시한다.
- 주요 파일: engines/native_gradient.py, engines/pyscf_engine.py(nuclear_gradient), engines/dispersion.py, molecular.py, main.py, molecular_app.py, workflow.py, 신규 test_gradient.py, scripts/validate_gradient.py, test_molecular_gui.py, README_KO.md, docs/NATIVE_DFT_KO.md·VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 검증: 고정 격자 유한차분 1.2e-9~3.3e-9 Ha/bohr(5조건 전 좌표), PySCF 기울기와 ≤8.8e-13, 격자 응답 크기 수준1~4 6.0e-4→1.1e-6, D4 기울기 FD 1e-9. test_gradient 6개, WSL 60개, Windows 38개+새 GUI 1개 통과. 상세 docs/VALIDATION_KO.md 0.8.3, outputs/gradient_validation/.
- 발견(미해결): 자체 DIIS가 근축퇴 열린 껍질 OH(거친 격자·엄격 수렴)에서 진동·미수렴. 기울기와 무관, 시험은 NH₂로 대체. 구조 최적화 전에 수정 필요.
- 한계: 주기계 힘/응력, 구조 최적화·제약·스캔, Hessian/진동은 미구현. ACTIVE_TASK.md 유지.

## QS-20261029-NATIVE-DFT · 단계 3 (DFT-D4 분산 보정) 완료 2026-10-29 · 0.8.2 · Claude

- 분자 계산에 `dispersion=d4`(GUI `분산 보정`, CLI `--dispersion d4`)를 추가했다. Grimme 그룹 공식 구현 dftd4로 D4 에너지(EEQ 전하 의존 C6, BJ 2체 + ATM 3체)를 계산해 SCF 후 총에너지에 별도 성분으로 더한다. 사용자 기준 "더 정확한 것"에 따라 D3 자체 구현 대신 D4 참조 구현을 택했다.
- 주요 파일: 신규 engines/dispersion.py, requirements-dispersion.txt, test_dispersion.py; molecular.py, engines/pyscf_engine.py, scf_control.py, energy_analysis.py, main.py, molecular_app.py, workflow.py(버전), test_basis_overrides.py(과거 서명 재현에 새 필드 제외), README_KO.md, START_HERE_KO.md, docs/NATIVE_DFT_KO.md·VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- PBE/PBE0/B3LYP만 허용(LDA 매개변수 없음), dftd4 지원 110원소(Rf~Rg 제외) 검증. 고립 분자는 라이브러리 기본 cutoff(2체 60 bohr)를 넓혀 모든 쌍 포함. 결과에 `scf_energy_hartree`, 성분 `dispersion`, 매개변수·EEQ 전하·cutoff·버전 기록. 체크포인트 서명에서 제외(밀도 무관). 반응 분석은 분산 모델·매개변수·버전 일치 요구. dftd4가 없어도 none 경로는 동작.
- 검증: test_dispersion 7개(2체 식 수동 재현 1e-10, dftd4 PySCF 어댑터와 1e-12, 불변성·전하·단위), WSL 전체 54개, Windows 37개 통과(7 skip 사유 기록). 물 PBE-D4 CLI 두 엔진 차 2e-12 Ha. 상세는 docs/VALIDATION_KO.md 0.8.2.
- 한계: D4 정확도 벤치마크를 재검증하지 않음, B3LYP VWN 변형 차이 가능성, 분산 힘·주기계 D4 미구현. ACTIVE_TASK.md 유지(힘/최적화 이후 단계 남음).

## QS-20261029-NATIVE-DFT · 단계 2 (원소별 기저) 완료 2026-10-29 · 0.8.1 · Claude

- 분자 계산에서 원소마다 다른 기저를 지정하는 `DFTConfig.basis_overrides`를 추가했다. GUI `원소별 기저` 입력란, CLI `--basis-element EL=BASIS`, 설정 JSON에서 사용한다. 기존 설정·결과·체크포인트와 호환된다.
- 주요 파일: molecular.py, engines/pyscf_engine.py, scf_control.py, energy_analysis.py, main.py, molecular_app.py, workflow.py(버전), 신규 test_basis_overrides.py, test_energy_analysis.py, test_molecular_gui.py, README_KO.md, docs/ELEMENTS_KO.md·NATIVE_DFT_KO.md·VALIDATION_KO.md·REQUIREMENTS_STATUS_KO.md.
- 원소별 실제 기저 데이터 검사, def2-ECP는 def2가 지정된 무거운 원소에만, 하나라도 Dyall이면 auto X2C/Gaussian/universal, ECP+X2C 혼합 차단. 결과에 `basis_by_element`, 혼합·미사용 지정 경고. 반응·흡착 분석은 같은 원소의 기저/ECP가 모든 항에서 같아야 한다(기존 `ecp_policy` 키 대체).
- 검증: 신규 10개(WSL) 통과, 기존 WSL 37개 통과, Windows GUI 실계산 포함 24개 통과(2개는 PySCF 필요로 WSL에서 통과). 혼합 기저 물 CLI 실계산 PySCF/자체 차이 3.1e-11 Ha. 상세는 docs/VALIDATION_KO.md 0.8.1.
- 한계: 혼합 기저 총에너지는 다른 조합과 직접 비교 불가, 정확도 검증 아님. 주기계 원소별 기저는 미지원. 분산 보정·힘/최적화 이후 항목은 미완료이며 ACTIVE_TASK.md를 유지한다.

## QS-20261029-NATIVE-DFT · 첫 구현 단계 완료 2026-10-29 · 0.8.0

- 첨부 프롬프트의 우선 구현 대상인 자체 분자 RKS/UKS SCF를 추가했다. 밀도·기울기·J/K·LDA/GGA XC 행렬·global hybrid·Pulay DIIS·혼합·총에너지 성분·엄격한 에너지/밀도/잔차 수렴을 직접 구현했다. 적분/AO/격자는 PySCF, XC 함수값·미분은 Libxc를 사용한다.
- 주요 파일: engines/native_scf.py, engines/pyscf_engine.py, molecular.py, scf_control.py, dft_workflow.py, dft_worker.py, molecular_app.py, main.py, geometry_tools.py, energy_analysis.py, workflow.py, examples/native_water.json, scripts/validate_native.py, test_native.py, test_energy_analysis.py, test_molecular_gui.py 및 안내/검증 문서.
- GUI/CLI 엔진 선택과 출처, 분자 중단·원자적 체크포인트·검사 후 재시작, 실시간 SCF 그래프, 구조 미리보기·거리/각도, 수치 프리셋, 실행 비교와 수지/조건을 검사하는 반응·흡착에너지 분석을 연결했다.
- 검증: WSL 물리·분석·기존 분자/주기계39개, Windows GUI10개(skip 없음), 기존 수소 workflow5개 모두 통과. 5종 실제 비교와 실제 CLI/에너지 분석 결과 보존. He/H2/H2O 독립 외부 대비 최대7.11e-14 Ha, OH 독립차2.40551e-6 Ha와 동일밀도 재시작차2.84e-14 Ha를 구분해 기록했다.
- 한계: 자체 H~Kr/비상대론/전체전자/점핵/160 AO, 외부의 기존118원소 경로 유지. 전체 원소·화합물의 정확도·바닥상태 또는 MLIP 대비 우위를 검증하지 않았다.
- **전체 첨부 요구사항 완료 기록이 아니다.** 후속 분산력·힘/최적화·제약/스캔·진동/열역학·NEB·응력/밴드/DOS 등은 미완료다. ACTIVE_TASK.md는 후속 인계를 위해 유지한다. 상세 상태와 증거는 docs/REQUIREMENTS_STATUS_KO.md 및 docs/VALIDATION_KO.md를 참조한다.

## QS-20261028-118 · 2026-10-28 완료 · 0.7.0

- 분자 입력을 H~Og 118개로 확장했다. 실제 기저 카탈로그, Dyall-v2z/v3z, scalar X2C, Gaussian 유한 핵, 원소별 def2 코어28/46/60, 초중원소 적분 격자를 추가했다.
- 주요 파일: molecular.py, element_catalog.py, data/element_capabilities.json, engines/pyscf_engine.py, molecular_app.py, molecular_plotting.py, main.py, periodic.py, test_molecular.py, test_molecular_gui.py, test_periodic.py, scripts/build_element_catalog.py, scripts/validate_heavy_elements.py, examples/auplus.json, examples/og_x2c.json. 사용·지원·검증 문서 및 환경 요구사항 갱신.
- GUI 원소별 지원표와 Au+/Og 예제, CLI 상대론·핵·격자 설정, X2C 출력 표기, 주기계 미지원 데이터 사전 차단을 추가했다.
- 검증: 분자·주기계30개 통과, 1,062개 원소/기저 조합 대조, Dyall 두 기저의 각118종 분자 구성 및118종 격자 생성, Windows GUI 입력4개 및 새 X2C 실계산1개, 기존 workflow5개 통과.
- Au⁺/Cs⁺/Ce⁴⁺/U⁶⁺/Og 실제 SCF·전자 수 검증. Og 방사 격자를 보강해 level2/3 에너지 차이6.09e-6 Ha 확인. 계산 조건·실패 원인·산출물 경로는 [검증 기록](docs/VALIDATION_KO.md), 지원 범위는 [원소 안내](docs/ELEMENTS_KO.md)에 보존.
- 한계: 118종 전체의 SCF·바닥상태·실험 정확도를 검증한 것은 아니다. 주기계 MOLOPT-SR은71원소, ML은 모델별 범위이며 SOC는 미구현이다.
- ACTIVE_TASK.md를 COMPLETED로 갱신하고 완료 기록·파일 보존을 확인한 뒤 임시 인계서를 삭제했다. 다음 작업에서는 새 인계서를 생성한다.

## QS-20261028-PROTOCOL · 2026-10-28 완료

- 사용자 지정 다중 AI 협업·작업 인계 규칙을 AGENTS.md에 등록하고 CLAUDE.md에서 연결했다.
- PROJECT_STATE.md에 시스템·환경·구현·연구 한계를, ACTIVE_TASK.md에 진행 중 원소 확장의 완료/미완료·실패 재현·재개 절차를 기록했다.
- 검증: 기존 관리 파일이 없음을 확인한 뒤 생성했으며, 실제 소스/결과 파일·테스트 출력·Git 상태와 대조했다. Git 저장소는 없으며 Git 이력 보존을 주장하지 않는다.
- 원소 확장은 별도 진행 작업 QS-20261028-118로 유지한다. 이번 규칙 등록 완료를 원소 확장 완료로 취급하지 않는다.
