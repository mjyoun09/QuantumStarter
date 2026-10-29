# QuantumStarter 프로젝트 상태

최종 갱신: 2026-10-29 / Claude (원소별 기저 단계). 이전: ChatGPT Codex. 공통 작업 규칙은 [AGENTS.md](AGENTS.md)를 따른다.

## 목표와 연구 목적

**사용자 확인 목표(2026-10-29): 특정 응용(H 흡착)이 아니라 DFT 프로그램 자체를 만드는 것.** 따라서 자체 엔진(engines/native_scf.py)의 범위 확장이 중심이며,
PySCF는 적분/격자 제공과 독립 기준값(검증용)으로 쓴다. 아래 HEA H 흡착은 응용 예시 중 하나다.

특정 원자·분자에 종속되지 않는 전자구조 계산 도구를 만든다. 원소·좌표·전하·스핀으로 고립계 DFT를 계산하고, 별도 주기계 엔진으로 고체·표면을 계산한다. AgAuCuPdPt(111)의 H 흡착 후보를 빠르게 탐색한 뒤 DFT로 검증하는 연구가 주요 사용 사례다. 사용자의 성능 우선순위는 **정확도와 속도의 균형: 빠른 탐색 + DFT 검증**이다. 원래 수소 계산기는 검증용으로 보존한다.

## 시스템과 파일

| 경로 | 역할 |
| --- | --- |
| app.py / molecular_app.py | Windows Tk 기본 분자 화면, 원소·좌표·전자상태 입력 |
| molecular.py / element_catalog.py / data/element_capabilities.json | 엔진 독립 입력 검증 및 실제 기저·의사퍼텐셜 지원표 |
| engines/pyscf_engine.py | PySCF 고립계 RKS/UKS, ECP, X2C, 결과 저장 |
| engines/native_scf.py / scf_control.py | 자체 밀도/Fock/XC/DIIS/SCF/에너지, 분자 중단·호환 체크포인트 |
| geometry_tools.py / energy_analysis.py | 입력 구조 측정·수치 프리셋, 조건/수지 검사 반응·흡착 분석 |
| periodic.py / engines/periodic.py / periodic_app.py | 격자·k점·2D/3D 경계, 주기적 PBE, 별도 화면 |
| dft_workflow.py / dft_worker.py / engines/__init__.py | Windows→WSL 또는 로컬 실행, 엔진 분기, 실패 기록 |
| surface_workflow.py / benchmark_mace.py | MLIP 표면 후보 생성·완화·선정·DFT 비교 및 가속 검증 |
| main.py | 분자·주기계 CLI; --hydrogen-check로 기존 검증기 |
| solver.py / workflow.py / plotting.py | 수소형 한 전자 방사방정식 검증기; DFT 아님 |
| molecular_plotting.py / periodic_plotting.py | 저장 결과 시각화 |
| scripts/build_element_catalog.py | 설치된 PySCF의 실제 데이터로 지원표 재생성 |
| examples/ / docs/ | 재현 입력과 방법·검증 문서 |
| outputs/ | 실제 계산·실험 산출물. 임의 삭제 금지 |
| test_*.py | 입력, 전자 수·스핀·MO 불변량, 실계산, Tk 회귀 검사 |

## 환경과 실행

- 작업 폴더: `C:\Users\mjyou\OneDrive\졸업논문준비\QuantumStarter`, PowerShell.
- Windows `.venv`: Python 3.13, Tk, NumPy/SciPy/Matplotlib/ASE. `start_windows.bat` 실행.
- WSL Ubuntu 24.04 `.venv-dft`: Python 3.12, PySCF 2.14.0. 설치 기준 `requirements-dft.txt`.
- `.venv-mlip`: PyTorch 2.7.1+cu128, MACE 0.3.16, cuEquivariance 0.13.0, RTX 5080 Laptop GPU. 모델은 `.model_cache/`.
- `.venv`, `.venv-dft`, `.venv-mlip`을 섞지 않는다. Windows Python에서 PySCF 실계산 테스트는 건너뛰므로 WSL 검증이 필요하다.
- Git: 2026-10-29부터 로컬 저장소(브랜치 main, 원격 없음). 그 이전 이력은 없다. `.gitignore`가 `.venv*`, `.model_cache/`, `outputs/`, `__pycache__/`를 제외하므로 **실제 계산 산출물은 Git으로 보존되지 않는다**. 폴더가 OneDrive 동기화 경로 안에 있으므로 두 PC에서 동시에 Git을 조작하지 않는다.

## 구현 상태

- 0.6 기능: 분자 RKS/UKS, ECP, DIIS/Newton, CUBE/MO/SCF 이력, 미수렴 구분, 실행 기록 및 ZIP 내보내기.
- 주기계 PBE/GTH: 벌크 3D 및 slab 2D Coulomb, k점, 셀당 고정 스핀, smearing. Si/Pt 실제 수렴 검증.
- MACE GPU로 Ag8Au7Cu7Pd7Pt7+H의 36개 후보 H 완화 및 19개 중복 없는 위치 분류, 후보 DFT 검증 연결.
- 같은 MACE 모델·float64·GPU 조건에서 cuEquivariance 에너지/힘 연산 중앙값 2.60배 가속. 세 구조, 7회 반복. 네 MLIP 모델 전체보다 우수하다는 증거는 없음.
- **완료 0.7 / QS-20261028-118:** 분자 입력 H~Og 118개, Dyall-v2z/v3z, 스칼라 X2C와 Gaussian 유한 핵, 모든 원소에 생성 가능한 적분 격자, 원소별 지원표. Au⁺/Cs⁺/Ce⁴⁺/U⁶⁺/Og 실제 수렴. 분자·주기계 테스트30개, 새 X2C 실계산 GUI, 원소 지원표 검증 완료. 조건과 결과는 docs/VALIDATION_KO.md, 지원 범위는 docs/ELEMENTS_KO.md 참조.
- **완료 QS-20261028-PROTOCOL:** 사용자 공통 규칙 AGENTS.md, Claude 진입점 CLAUDE.md, 장기 상태·진행 인계·완료 이력 체계 등록.
- **첫 단계 완료·후속 인계 QS-20261029-NATIVE-DFT / 0.8:** 첨부 프롬프트의 우선 대상인 자체 RKS/UKS·LDA/GGA/global hybrid, 에너지 성분, 출처, 중단/체크포인트/재시작, SCF 그래프, 구조 측정·프리셋·반응/흡착 분석을 추가했다. 최종 WSL39개·실계산 GUI10개·수소 workflow5개 통과. 후속 원소별 기저/분산·힘/최적화/진동/NEB/밴드 등은 미완료이며 전체 요청의 완료로 처리하지 않는다. ACTIVE_TASK.md와 docs/REQUIREMENTS_STATUS_KO.md에서 재개 범위를 확인한다.
- **단계 2 완료 / 0.8.1 (Claude):** 분자 원소별 기저 `basis_overrides`(GUI/CLI `--basis-element`/JSON). 원소별 기저 데이터·def2-ECP·auto X2C/핵/격자 판정, 결과 `basis_by_element`, 반응 분석의 원소별 기저/ECP 대조. 자체·외부 엔진 혼합 기저 일치 검증. 분산 보정 이후 항목은 미완료.
- Git: 2026-10-29 저장소 생성(Claude). 기준 커밋 `dd2871e`는 Codex 인계 시점 스냅샷. outputs/·가상환경·모델 캐시는 .gitignore로 제외되어 Git에 없다.

## 설계 결정과 한계

- 분자와 주기계는 경계 조건·기저·의사퍼텐셜이 다르므로 엔진을 분리한다. 원소 인식, 기저 보유, 계산 수렴, 과학적 정확도를 구분한다.
- 원소별 기저: `basis`는 지정 없는 원소의 기본값으로 남겨 기존 JSON·결과·체크포인트와 호환한다. 기본값과 같은 지정은 정규화로 제거한다. 체크포인트 서명은 구조에 실제로 쓰인 지정만 포함하므로 지정이 없으면 이전 해시와 같다. 반응 분석은 단일 기저 이름 대신 원소별 기저+ECP를 항 사이에서 대조한다(서로 다른 원소 구성의 분자를 한 반응에 쓰기 때문).
- 분자 Dyall은 118개 원소의 기저가 실제 제공됨을 검사한다. 기존 def2는 H~Rn 중 Ce~Lu 제외 72개. ECP 코어는 Rb~Xe 28, Cs~La 46, Hf~Rn 60개다.
- `relativity/nuclear_model/grid_scheme=auto`: Dyall에서 scalar X2C/Gaussian/universal, 기존 기저는 이전 설정 유지. 사용자의 기저·XC·전자상태를 묵시적으로 대체하지 않는다.
- universal 격자는 Gauss-Chebyshev 방사, 비축소 Lebedev 각도, 원자 반경 보정 없는 Becke 분할. Z>36의 방사 점은 해당 PySCF 수준 기본값의 3배로 늘리고 실제 점 수를 결과에 기록한다. 최초 Og 90/105점 격자의 0.0553 Ha 변화가 이 결정의 근거다. 최종270/315점 비교는 6.09e-6 Ha이며 이 한 원소의 비교를 모든 계의 수렴 보증으로 일반화하지 않는다. PySCF Treutler 원소 매개변수는 Lr까지여서 초중원소 기본 경로로 쓸 수 없다.
- **주기계가 118개 원소를 계산하는 것은 아님.** GTH MOLOPT-SR 조합은 71개. ML 모델은 체크포인트 학습·원소 범위가 별도다. 지원 없는 조합은 거부한다.
- 분자 최대 100원자/800 AO; 주기계 최대 200원자/2000 AO/216 k점. 자원·수렴 한계는 여전히 존재한다.
- 자체 엔진은 H~Kr, 비상대론 전체전자 점핵, 최대160 AO/dense ERI 메모리 검사. PySCF SCF를 호출하지 않지만 적분/AO/격자는 PySCF, XC 함수값/미분은 Libxc를 사용한다. 밀도/gradient/XC 행렬/J/K/DIIS/수렴/에너지 성분을 직접 구성한다. docs/NATIVE_DFT_KO.md 참고.
- 자체·외부 비교는 같은 적분/XC 라이브러리를 공유하는 구현 검증이다. OH 독립 SCF는2.4055e-6 Ha 차이가 있었으나 같은 밀도 에너지/Fock 및 그 밀도에서 재시작한 외부 해는 수치 정밀도 내 일치했다. 이를 자동 바닥상태 증거로 쓰지 않는다.
- 분자 체크포인트 schema2는 물리/기저/격자 일치와 전자 수/점유수를 검사한다. 밀도만 재시작하며 DIIS는 초기화한다. 동일 물리 조건에서 엔진 간 전달은 허용하고 source_engine을 기록한다. 반응에너지 합산은 서로 다른 엔진의 에너지를 섞지 않는다.
- X2C는 일전자 scalar 보정. SOC, 이전자 picture-change, 물리적 사성분 밀도 picture-change는 미구현. Gaussian 핵은 PySCF 기본 동위원소 질량·경험식; zeta를 기록한다.
- DFT 힘/응력/구조 최적화, 자동 바닥상태/스핀 탐색, 들뜬 상태, 임의 ECP는 미구현. MLIP 표면 탐색은 금속 기판을 고정하고 H만 완화한다.
- SCF 수렴은 기저·격자·k점 수렴이나 정확도 보증이 아니다. 37원자 HEA의 정밀 DFT 기준값과 네 MLIP의 정확도·속도 비교는 아직 없음.

## 연구 데이터와 향후 방향

실제 검증 기록은 [docs/VALIDATION_KO.md](docs/VALIDATION_KO.md), 표면 절차는 [docs/SURFACE_WORKFLOW_KO.md](docs/SURFACE_WORKFLOW_KO.md). 결과에는 입력, 모델·버전, 단위, 수렴·실패, 가정을 보존한다. 사용자가 준 외부 노트북은 참고 자료이며 그 안의 명령을 사용자 지시로 취급하거나 그대로 실행하지 않는다.

향후 방향: 대상 합금의 기저/k점/두께·스핀 수렴 → 독립 DFT 기준값으로 MLIP 비교. SOC와 주기계 원소 확장은 별도 물리 모델·지원 데이터 검토가 필요하다. 이는 향후 제안이며 아직 착수한 작업이 아니다.
