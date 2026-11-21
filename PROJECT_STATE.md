# QuantumStarter 프로젝트 상태

최종 갱신 2026-11-21 / ChatGPT Codex. 원자 SCF에 Anderson 혼합을 추가했다(0.9.8, QS-20261121-ATOM-MIX). 기본 선형 혼합과 비교해 H–Kr PBE 72조건 모두 기존 기준으로 수렴했고 반복 수 합계가 3,768회에서 1,060회로 줄었다. [검증 결과](docs/VALIDATION_KO.md) · [공통 작업 규칙](AGENTS.md)

같은 날 QS-20261121-N2-THERMO에서 중성 질소의 평형 구조·조화 진동수로 1 bar, 298.15–2,000 K의 열물성을 계산해 NIST와 대조했다. 열용량은 36개 온도에서 모두 기준보다 작았고 최대 편차는 2,000 K에서 −0.292 J/(mol·K)였다. [질소 검증 기록](docs/VALIDATION_KO.md) · [온도별 표](docs/validation/n2_thermo.csv)

QS-20261121-H2-DISS에서 고정 거리 H₂의 PBE/def2-SVP 해리 곡선을 검증했다. 19개 거리의 RKS·UKS 에너지가 독립 PySCF SCF와 일치했고 스핀 곡률이 1.6178125–1.618125 Å에서 부호를 바꿨다. 이 H₂ 연구 결과는 장기 미구현 항목인 일반 분자 RKS→UKS 안정성 진단을 개발할 때 쓸 기준이다. [계산 조건과 범위](docs/VALIDATION_KO.md) · [거리별 표](docs/validation/h2_scan.csv)

## 목표와 연구 목적

**사용자 확인 목표(2026-10-29): 특정 응용(H 흡착)이 아니라 DFT 프로그램 자체를 만드는 것.** 따라서 자체 엔진(engines/native_scf.py)의 범위 확장이 중심이며,
PySCF는 적분/격자 제공과 독립 기준값(검증용)으로 쓴다. 아래 HEA H 흡착은 응용 예시 중 하나다.

**사용자 확인 응용 방향(2026-10-30): HEA 표면 응용이 목적이 아니다. DFT를 항공우주와 엮어 물리 쪽으로 활용하고 싶다.** 구체 주제는 아직 미정(지도교수 확인 필요). 후보와 현재 엔진 적합도는 대화에서 제시했다: 고온 기체·추진제 화학종의 통계역학 열물성(Cp/H/S, NASA 다항식)과 반응 속도(TST)가 현재 엔진(경원소 분자, TS/NEB, 진동수·RRHO)에 가장 잘 맞는다. 이후 HEA 관련 기록은 과거 응용 예시로만 본다.

특정 원자·분자에 종속되지 않는 전자구조 계산 도구를 만든다. 원소·좌표·전하·스핀으로 고립계 DFT를 계산하고, 별도 주기계 엔진으로 고체·표면을 계산한다. AgAuCuPdPt(111)의 H 흡착 후보 탐색·DFT 검증은 초기 응용 예시이며 현재 연구 방향은 위 항공우주 물리다. 사용자의 성능 우선순위는 **정확도와 속도의 균형: 빠른 탐색 + DFT 검증**이다. 원래 수소 계산기는 검증용으로 보존한다.

## 시스템과 파일

| 경로 | 역할 |
| --- | --- |
| app.py / molecular_app.py | Windows Tk 기본 분자 화면, 원소·좌표·전자상태 입력 |
| molecular.py / element_catalog.py / data/element_capabilities.json | 엔진 독립 입력 검증 및 실제 기저·의사퍼텐셜 지원표 |
| engines/pyscf_engine.py | PySCF 고립계 RKS/UKS, ECP, X2C, 결과 저장 |
| engines/native_scf.py / scf_control.py | 자체 밀도/Fock/XC/DIIS/SCF/에너지, 분자 중단·호환 체크포인트 |
| engines/native_newton.py / native_stability.py / native_response.py | 자체 Newton·실수 내부 안정성, 해석적 J/K/XC 선형응답 |
| engines/optimizer.py / irc.py / neb.py / frequencies.py | geomeTRIC 최적화·제약·스캔·TS, IRC(끝점 판정·최소화·경로 복구), CI-NEB, 유한차분 진동수·RRHO |
| geometry_tools.py / energy_analysis.py | 입력 구조 측정·수치 프리셋, 조건/수지 검사 반응·흡착 분석 |
| periodic.py / engines/periodic.py / periodic_app.py | 격자·k점·2D/3D 경계, 주기적 PBE, 별도 화면 |
| engines/periodic_optimizer.py / periodic_convergence.py | FFTDF bulk 원자·셀 완화, 고정 구조 FFT mesh 수렴 검사 |
| dft_workflow.py / dft_worker.py / engines/__init__.py | Windows→WSL 또는 로컬 실행, 엔진 분기, 실패 기록 |
| surface_workflow.py / benchmark_mace.py | MLIP 표면 후보 생성·완화·선정·DFT 비교 및 가속 검증 |
| main.py | 분자·주기계 CLI; --hydrogen-check로 기존 검증기; --model he-lda로 방사 격자 He LDA; --model atom으로 H–Kr 원자 LDA/LSDA 및 PBE |
| solver.py / workflow.py / plotting.py | 수소형 한 전자 방사방정식 검증기; DFT 아님. solver.radial_grid·lowest_radial_states는 atomic_dft와 공유 |
| atomic_dft.py / atomic_plotting.py | 중성 He 구대칭 Kohn–Sham LDA-SCF(균일 방사 격자, Slater+VWN5 닫힌 식, 하트리·SCF·에너지 자체 구현, PySCF·Libxc 런타임 미사용). atomic_plotting은 원자 그림(create_atom_figure)도 담당 |
| radial_atom.py | H–Kr 구대칭 원자 LDA/LSDA(로그 격자 r = a(eˣ−1), Liouville 변환 대칭 삼중대각, 여러 껍질·사용자 배치, VWN5 스핀 보간 자체 구현) |
| atomic_pbe.py / scripts/validate_atomic_pbe.py / scripts/validate_spherical_gaussians.py | 자체 PBE와 이산 변분 퍼텐셜, PySCF 및 독립 구대칭 가우시안 기준 계산 |
| data/nist_atomic_reference.json / scripts/build_nist_atomic_reference.py | NIST Atomic Reference Data(Z=1–36, LDA·LSD) 원 페이지 파싱 결과와 생성기(페이지 SHA-256 포함) |
| molecular_plotting.py / periodic_plotting.py | 저장 결과 시각화 |
| scripts/build_element_catalog.py | 설치된 PySCF의 실제 데이터로 지원표 재생성 |
| examples/ / docs/ | 재현 입력과 방법·검증 문서 |
| outputs/ | 실제 계산·실험 산출물. 임의 삭제 금지 |
| test_*.py | 입력, 전자 수·스핀·MO 불변량, 실계산, Tk 회귀 검사 |

## 환경과 실행

- 작업 폴더: `C:\Users\mjyou\OneDrive\졸업논문준비\QuantumStarter`, PowerShell.
- Windows `.venv`: Python 3.13, Tk, NumPy/SciPy/Matplotlib/ASE. `start_windows.bat` 실행.
- WSL Ubuntu 24.04 `.venv-dft`: Python 3.12, PySCF 2.14.0. 설치 기준 `requirements-dft.txt`. 선택: dftd4 4.3.1 + tomli 2.5.0 (`requirements-dispersion.txt`), geomeTRIC 1.1.1 + networkx 3.7 (`requirements-optimize.txt`); 둘 다 2026-10-29 이 PC에 설치됨.
- `.venv-mlip`: PyTorch 2.7.1+cu128, MACE 0.3.16, cuEquivariance 0.13.0, RTX 5080 Laptop GPU. 모델은 `.model_cache/`.
- `.venv`, `.venv-dft`, `.venv-mlip`을 섞지 않는다. Windows Python에서 PySCF 실계산 테스트는 건너뛰므로 WSL 검증이 필요하다.
- **사용자 요청(2026-10-29): 검증 중 사용자 화면에 창을 띄우지 않는다.** 분자·주기계 GUI 시험은 `QUANTUMSTARTER_TK_HIDDEN=1`(루트·Toplevel을 매핑 전에 숨김, 화면 배치 검사만 생략)로 실행한다. 주기계 숨김 모드는0.9.3에 추가됐고 실계산은 `QS_PERIODIC_GUI_DFT=1`로 활성화한다. 숨김 모드가 없는 test_gui(수소 화면)는 사용자 요청 없이 실행하지 않는다. 분자 예: `$env:QUANTUMSTARTER_TK_HIDDEN='1'; $env:QUANTUMSTARTER_GUI_DFT_TEST='1'; .\.venv\Scripts\python.exe -m unittest test_molecular_gui`.
- **줄바꿈 주의**: 저장소 파일은 LF다. Windows `.venv` Python의 `Path.write_text`는 CRLF로 바꿔 쓰므로 스크립트로 파일을 고칠 때는 bytes로 쓰거나 `newline='\n'`을 지정한다. Git Bash의 `grep -c $'\r$'`은 CRLF 검사로 믿을 수 없다(바이트 검사 사용). 예외: examples/native_water.json은 처음부터 CRLF로 커밋됨.
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
- **단계 3 완료 / 0.8.2 (Claude):** 분자 DFT-D4 분산 에너지(`dispersion=d4`), 공식 dftd4 라이브러리(선택 의존성, requirements-dispersion.txt). PBE/PBE0/B3LYP, 110원소, 고립 분자는 절단 없이 합산. 결과에 SCF/분산 분리·매개변수·EEQ 전하 기록, 반응 분석은 분산 조건 일치 검사. 분산 힘·주기계 D4 미구현. 사용자 선택 기준: "더 정확한 것" → 모델이 더 정확한 D4, 참조 구현 사용.
- **단계 4 완료 / 0.8.3 (Claude):** 분자 원자 힘 `calculation=gradient`. 자체 엔진은 engines/native_gradient.py의 해석적 RKS/UKS 기울기(LDA/GGA/global hybrid, 성분별 저장), PySCF 엔진은 nuc_grad_method(ECP/X2C/Newton), D4 기울기 합산. 고정 격자 유한차분 ~1e-9, PySCF와 ~1e-12 Ha/bohr 일치. 격자 응답 항 미포함(수준 3에서 ~1e-5 Ha/bohr).
- **단계 5 완료 / 0.8.4 (Claude):** 자체 DIIS 정체 수정(de7d245)과 `calculation=optimize`(geomeTRIC TRIC + 자체/PySCF 해석적 힘, 단계별 기록, 이전 구조 밀도 초기 추정). PySCF 자체 geomeTRIC 드라이버와 최적 구조 소수 6자리 일치. 사용자 선택: 최적화기 geomeTRIC.
- **단계 6 완료 / 0.8.5 (Claude):** `calculation=frequencies` — 해석적 힘 중심 유한차분 Hessian, 정확한 동위원소 질량(COMMON_ISOTOPE_MASSES; PySCF isotope_avg=False는 정수 질량수라 사용 안 함), 병진·회전 정확 제거, 허수 진동 판정, ZPE·진동 열역학. PySCF 해석적 Hessian과 ≤0.5 cm⁻¹. 함께 수정: 실행 폴더 이름에 마이크로초(같은 초 실행의 기록 순서 버그).
- **단계 7-B 완료 / 0.8.6 (Claude):** 자체 SCF 원자 밀도 중첩 초기 추정·잔차 기준 혼합(반복 30~50% 감소, Fe 수렴), 자체 SCF 안정성 분석(실수 내부, check/follow, PySCF 판정과 일치). 사용자 지시 "전부 해" → 이어서 7-A(자유에너지), 7-C(제약·스캔·NEB/TS).
- **단계 7-A 완료 / 0.8.7 (Claude):** RRHO 이상기체 열역학(ZPE/H/G, 대칭수 자동·입력, 다중도), 반응 분석의 ZPE/H/G 선택과 온도·압력 일치 검사. PySCF thermo와 H·G 1e-10 상대 일치.
- **단계 7-C 완료 / 0.8.8 (Claude):** geomeTRIC 제약 최적화, 1차원 이완 스캔, P-RFO 전이상태 탐색(진동수로 1개 허수 확인). NEB·IRC는 미착수(geomeTRIC 1.1.1 API 존재 확인).
- **단계 8 완료 / 0.8.9 (Claude, 2026-10-30):** CI-NEB(geomeTRIC 사슬 최적화 + 자체/PySCF 힘, 선형 보간, 양 끝 고정). NH₃ 반전에서 TS 탐색 장벽과 0.001 kJ/mol 일치. 사용자 지시 "전부 해"(B·A·C·NEB) 완료. IRC는 미구현.
- **단계 9 / 0.9.0 (Claude, 2026-10-30, 사용자 "둘 다 해볼래"):** 자체 신뢰영역 Newton SCF(안장점 탈출 포함; 물·Fe·OH 검증, Cr·Co는 여전히 미수렴), 주기계 DOS·밴드(Si 검증; 힘·응력은 PySCF GDF 미지원으로 보류).
- **단계 10 / 0.9.1 (ChatGPT Codex, 2026-10-30):** 자체 해석적 궤도 Hessian 응답(Libxc f_xc + 자체 J/K/XC 조립)·궤도 간격 전처리·CG 정밀도 및 신뢰영역 수락 개선. 기록된 Cr(스핀6)·Co(스핀3) PBE/STO-3G/grid3 정수 점유 사례 수렴 및 실수 내부 안정성 확인, 같은 밀도에서 독립 PySCF Newton 에너지 일치. 물 응답 7회 중앙값 FD 대비3.74배(전체 SCF 속도 보장 아님). graphene slab 밴드 실계산 검증. bulk FFTDF 힘·응력은 별도 검증 스크립트에서 에너지 FD 일치 확인했으며 제품 API 통합은 미구현. 상세 조건·로그는 docs/VALIDATION_KO.md와 ACTIVE_TASK.md.
- **단계 11 / 0.9.2 (Claude, 2026-10-30):** bulk FFTDF 주기계 힘·응력 제품 API(`integral_scheme=fftdf`, `calculation=gradient`). Γ·일반 k점·UKS·smearing(F 미분)을 유한차분으로 확인. slab은 GDF만이라 힘·응력 없음. 셀 최적화는 다음 후보.
- **단계 12 / 0.9.3 (ChatGPT Codex, 2026-10-30):** `cell_optimize` — bulk FFTDF 원자·셀 동시 완화(ASE BFGS/FrechetCellFilter), full/isotropic, 고정 FFT 정수 mesh·상대 k점, E 또는 smearing F 최소화. 힘·응력·ΔE/F 개별 판정과 최종 독립 SCF 재검사. FFT mesh 수렴 도구(연속 두 구간 에너지·힘·응력 비교), GUI/CLI와 주기계 취소 지원. Si 2³ k점 예제 등방7/전단 포함 full10단계 수렴, 초기·최종 mesh 검사 통과. 물성 정확도 검증은 아님. 사용법: docs/PERIODIC_CELL_KO.md.
- **단계 13 / 0.9.4 (Claude, 2026-11-01, 사용자 "해볼래"로 직전 제안 승인):** 분자 `calculation=irc` — geomeTRIC 1.1.1 내장 IRC + 자체/PySCF 힘. 끝점별 기울기 판정, 기준 미달 끝점은 일반 최소화로 마무리, 마지막 방향 단계 한도 시 geometric_irc.xyz로 경로 복구. NH₃ 반전(거울상 끝점)·HCN⇌HNC(끝점 = 독립 최적화 이성질체, 진동수 최소점) 검증. RRHO 열역학을 ASE IdealGasThermo와 교차검증(H·S ≤4.1e-7 상대, 298~3000 K). 상세는 docs/VALIDATION_KO.md 0.9.4.
- **QS-20261112-HE-LDA / 0.9.5 (Claude, 2026-11-12, 사용자 지시 "중성 He LDA-SCF 직접 구현·검증"):** atomic_dft.py 방사 격자 LDA-SCF(`main.py --model he-lda`). 기본(dr 0.001, r_max 30) E = −2.834832195 Ha, ε₁ₛ = −0.570423731 Ha, PySCF RKS 기준(−2.834835624, 짝수 템퍼드 s 39개)과 +3.4e-6 Ha. 격자 2차 수렴, Richardson은 PySCF와 2e-11 Ha, NIST LDA 표(−2.834836)와 6자리 일치. 상세 docs/VALIDATION_KO.md 0.9.5.
- **QS-20261117-ATOM-LSDA / 0.9.6 (Claude, 2026-11-17, 사용자 "진행해" — 직전 제안 "방사 격자 엔진을 임의 원자로 확장" 승인):** radial_atom.py(`main.py --model atom --element X [--spin-polarized]`). H–Kr 72건(LDA·LSD) 기본 격자(h 0.001) 모두 수렴; NIST 대비 기본 H–Ar |ΔE| ≤ 2.4e-5, K–Kr 상대 ≤ 5.7e-8; Richardson은 총에너지·성분·고윳값 모두 ≤ 8e-7(NIST 6자리 반올림 폭). 닫힌 껍질 5개 PySCF RKS와 ≤ 4.3e-8. 상세 docs/VALIDATION_KO.md 0.9.6. 다음 후보: GGA(PBE) 원자, 그다음 비구면·다중항 또는 의사퍼텐셜.
- **QS-20261117-ATOM-GGA / 0.9.7 (ChatGPT Codex, 2026-11-18)**. `--functional pbe`로 H–Kr 72건 모두 수렴했다. 독립 에너지 대조는 닫힌 껍질 7종과 열린 껍질 B에서 1e−6 Ha 기준을 통과했다. 기존 LDA/LSDA 72건은 같은 환경에서 원본과 비트 단위로 일치한다. [방법과 범위](docs/ATOMIC_PBE_KO.md)
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
- 분자 DFT 원자 힘(격자 응답 항 없음), geomeTRIC 구조 최적화, 유한차분 핵 Hessian·조화 진동수·RRHO, 실수 내부 SCF 안정성은 구현. 해석적 **핵 좌표** Hessian(궤도 Hessian 응답과 다름), quasi-RRHO·용매 자유에너지, 주기계 IRC, 자동 바닥상태/스핀 탐색, 삼중항·복소 SCF 안정성, 들뜬 상태, 임의 ECP는 미구현. 주기계 FFTDF bulk 힘·응력·셀 완화는 제품 API/GUI로 구현했으며 GDF/2D 경로에서는 미지원이다. PySCF 2.14 응력은 일부 적분 내부 유한차분을 포함한다. MLIP 표면 탐색은 금속 기판을 고정하고 H만 완화한다.
- bulk 셀 최적화는 외부압력0, 원자 제약 없음, 초기 정수 FFT mesh 고정. smearing F를 미분하는 힘·응력에 맞춰 F를 최소화한다. 등방 완화는 평균 응력만 판정한다. 중간 SCF 실패를 거부하며 최종 재계산은 `final/`에 두어 루트에 잘못된 수렴 상태를 남기지 않는다. 요청 설정과 최종 구조 설정을 구분한다. 최적화 재시작·포논 안정성·자동 k점/기저 수렴은 미지원.
- 2026-10-29 발견·수정: 자체 DIIS 계수 상한 30이 수렴 정체를 일으킴 → 의사역행렬로 수정(0.8.4). 근축퇴 열린 껍질(OH)은 다중 UKS 해가 있고 수렴한 해가 불안정(안장점)일 수 있음 → 0.8.6 자체 안정성 분석(check/follow)으로 검출·추적 가능. Fe 원자는 0.8.6에서 수렴, 기록된 **Cr·Co 정수 점유 미수렴은 0.9.1 Newton에서 해결**. 모든 기저/스핀/격자 및 전이금속 화합물의 수렴·바닥상태 보증은 아니다.
- IRC(0.9.4) 구동기는 Sella(미설치, jax 추가 의존) 대신 이미 쓰는 geomeTRIC 내장 IRC를 택했다. geomeTRIC은 한 방향이 단계 한도에 걸려도 전체 성공을 보고하고(HCN⇌HNC의 선형 HNC 끝에서 실측), customengine 분자에는 IRC 전 토폴로지가 없어 실패하므로, 성공 판정은 우리 끝점별 기울기 검사로 하고 토폴로지는 미리 만든다. 경로 s의 질량은 geomeTRIC 표준 원자량(진동 해석의 최빈 동위원소와 다름). IRC 끝점은 연결의 필요조건이며 의도한 반응물·생성물과의 대조는 사용자 몫이다.
- 방사 격자 He LDA(0.9.5): 분자 엔진(가우스 기저, PySCF 적분)과 별개인 실공간 경로로, Windows `.venv`만으로 실행된다. XC는 Libxc와 같은 정의의 닫힌 식(Libxc 7.0.0과 1e-14 일치를 시험)이라 런타임 의존성이 없다. 1s 고윳값은 고유벡터의 Rayleigh 몫으로 계산한다: eigh_tridiagonal 고유값의 절대 정밀도가 eps·2/dr²(dr 0.001에서 4e-10 Ha)여서 SCF ΔE가 요동했기 때문이다. 같은 이유로 밀도 잔차 하한이 ~0.4·eps·2/dr²이며 grid_points 상한은 1,000,000. 이 경로는 중성 He 1s² 전용으로 0.9.5 검증 기록 보존을 위해 그대로 둔다.
- 원자 LDA/LSDA(0.9.6, radial_atom.py): 로그 격자의 대칭 삼중대각 C는 원점 근처 성분 ~1/(h·a)²의 계단형이라 eigh_tridiagonal 기본 허용오차로는 고윳값이 Hartree 단위로 틀린다(실측) → stebz에 tol = 4·tiny를 주어 상대 정밀도로 풀고 ε는 Rayleigh 몫. 상대 이산화 오차 ≈ 0.05·h²(Z에 거의 무관)이므로 기본 h = 0.001에서 총에너지 상대 오차 ~5e-8. 구면 평균·Hund 점유는 NIST 표와 같은 관례(N LSD 고윳값 6자리 일치로 확인). Libxc 7은 1e-15 bohr⁻³ 미만 스핀 밀도를 1e-15로 바꿔 계산하므로 Cr LSD 꼬리에서 소수 스핀 퍼텐셜이 Libxc와 다르다(에너지 영향 ~1e-15 Ha, 자체 값은 한쪽 차분과 일치). 0.9.7은 같은 비상대론·점핵·구면 평균 모델에 PBE를 추가했다. 상대론·다중항·음이온은 미구현이다.
- SCF 수렴은 기저·격자·k점 수렴이나 정확도 보증이 아니다. 37원자 HEA의 정밀 DFT 기준값과 네 MLIP의 정확도·속도 비교는 아직 없음.

## 연구 데이터와 향후 방향

실제 검증 기록은 [docs/VALIDATION_KO.md](docs/VALIDATION_KO.md), 표면 절차는 [docs/SURFACE_WORKFLOW_KO.md](docs/SURFACE_WORKFLOW_KO.md). 결과에는 입력, 모델·버전, 단위, 수렴·실패, 가정을 보존한다. 사용자가 준 외부 노트북은 참고 자료이며 그 안의 명령을 사용자 지시로 취급하거나 그대로 실행하지 않는다.

향후 방향은 위 사용자 연구 목적의 **DFT × 항공우주 물리**를 우선한다. N₂의 온도별 Cp/ΔH/S 및 기준 대조를 첫 사례로 수행했다. 다른 기체로의 확대, NASA 7계수 적합, TST 속도는 후속 후보이며 구체 논문 주제는 확정 전이다. 코드 후속 후보는 quasi-RRHO·용매, 삼중항·복소 SCF 안정성이다(IRC는 0.9.4 완료). 원자 PBE는 0.9.7, Anderson 혼합은 0.9.8에서 완료했다. 원자 후속 후보는 비구면·다중항, 스칼라 상대론, 의사퍼텐셜이며 우선순위는 미정이다. 기존 합금/MLIP 예제는 보존하지만 사용자 최종 연구 대상으로 간주하지 않는다. SOC와 주기계 원소 확장은 별도 물리 모델·지원 데이터 검토가 필요하다.
