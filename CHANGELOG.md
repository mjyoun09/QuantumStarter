# 변경 이력

완료된 작업만 기록한다. 진행 작업은 ACTIVE_TASK.md, 장기 상태는 PROJECT_STATE.md를 참조한다.
규칙 등록 이전의 검증 이력은 [docs/VALIDATION_KO.md](docs/VALIDATION_KO.md)에 보존되어 있다.

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
