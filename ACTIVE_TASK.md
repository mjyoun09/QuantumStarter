# 현재 작업

## QS-20261124-MIX-REFERENCE · 로드맵 03 혼합기체 비열·음속의 독립 기준 대조

**작업 식별 정보** · ChatGPT Codex · COMPLETED · 2026-11-24 · 0.12.1 · 시작 HEAD `fe13aff`

**작업 목적** · 몰분율 N₂ 0.79/O₂ 0.21, 1 bar, 298.15–2,000 K에서 자체 DFT 유래 NASA7의 Cp와 음속을 독립 종별 열화학 자료로 산출한 값과 대조한다. 표에 기준의 종류·단위·온도 범위·최대 편차 및 그 온도를 남긴다. 원문에서 조건이 일치하는 직접 음속 측정점도 확인되면 별도로 대조한다.

**작업 계획**

1. NIST 원문·기존 Shomate 계수와 구간을 대조하고, 정확한 조성·압력·온도와 음향 조건을 확인할 수 있는 측정 자료를 조사한다.
2. `gas_mixture.py`의 기존 계산/내보내기와 `molecular_studies_app.py`의 혼합기체 화면에 기준 Cp·유도 음속과 편차·극값을 연결한다. CSV·그림·출처와 재현 절차를 보존한다.
3. 순수 기체·혼합 조성, 기준 온도 구간 경계, 독립 수치와 기본 36점 전 구간, CLI와 숨김 Windows GUI를 검증한 뒤 로드맵·장기 상태·변경 이력을 갱신한다.

**진행 상황** · 03의 GUI·CLI, 36점 기준표, 최대 편차표와 그림을 완성했다. 순수 기체·삼원 조성·온도 구간 경계, JANAF 표, 사용자 계수 교체, GUI 재계산을 검증했다. 후속 04–05는 대기 중이다.

**구현 및 실행 기록**

- `gas_mixture.py`가 NIST Shomate 기준 Cp와 유도 음속을 계산한다. `molecular_studies_app.py`와 CLI에서 부호 있는 편차를 확인할 수 있다. 기본 2,000 K의 편차는 −16.353966 J/(kg·K), +1.690953 m/s다. 검토한 이원계 실험은 범위보다 저온이었고 정밀 실험의 합성 공기에는 Ar이 포함되었다. [수식·출처·측정 조건](docs/MIXTURE_REFERENCE_KO.md)
- Windows 기본 테스트 231건 중 139건 통과·92건 건너뜀. 관련 테스트 15건은 Windows·WSL에서 각각 모두 통과했고, 별도의 Windows 숨김 GUI 3건도 통과했다. Cantera 3.2.0으로 같은 NASA7 계수의 혼합 연산 190조건을 재검사했다.
- `main.py gas-mixture --output outputs/mixture_reference_validation/release` 실행 결과의 두 CSV와 PNG는 `docs/validation/air_reference_*` 세 파일과 바이트 단위로 일치한다. 이원계·삼원계 실제 CLI 결과는 `outputs/mixture_reference_validation/`에 보존했다. [온도별 대조표](docs/validation/air_reference_comparison.csv) · [극값·적용 범위](docs/validation/air_reference_extrema.csv)

**다음 AI를 위한 인계** · [로드맵 04](ROADMAP.md)부터 시작한다. 먼저 `scripts/validate_air_sensitivity.py`와 `docs/validation/air_basis_grid.csv`에서 기존 기저×격자 계산 조건을 확인한다. 03의 NIST 기준을 고정하고 헤시안 변위·교환상관 함수 등 한 조건씩 바꾼 비교쌍에서 종별 진동수·Cp와 혼합 음속을 추적한다. 아래 별도 장기 HANDOFF_READY 인계가 있어 이 파일을 보존한다.

## QS-20261123-INFRARED · 로드맵 02 분자 적외선 세기와 동위원소 스펙트럼

**작업 식별 정보** · ChatGPT Codex · COMPLETED · 2026-11-23 · 0.12.0 · 시작 HEAD `0aa2caa`

**작업 목적** · 물과 중수의 조화 진동수·적분 IR 세기(km/mol)를 모드별로 계산하고 같은 축에서 비교한다. 완료 기준은 CLI·Windows GUI와 그림·수치 저장, H₂ 비활성 및 물 활성 모드, 변위 간격 민감도와 독립 계산 검증, 구버전 결과 호환이다.

**작업 계획**

1. `engines/pyscf_engine.py`·`engines/frequencies.py`에서 기존 ±변위 SCF에 사용한 총 밀도로 고정 원점 쌍극자를 구하고 중앙차분 미분을 저장한다. 지원하는 전자구조 모델과 단위를 명시한다.
2. `molecular_isotopes.py`와 기존 비교 탭·CLI에서 질량 정규화한 모드에 미분을 투영한다. 각 피크의 적분 세기, 가우스 선폭, 겹친 스펙트럼을 저장한다. 구버전 헤시안으로는 기존 진동수·ZPE 비교를 계속 제공한다.
3. 자체 DFT 물·H₂ 실계산, 독립 구현·변위 간격 대조, 숨김 GUI를 검증한다. 재현 가능한 물 예제와 검증 수치, 로드맵 및 변경 기록을 갱신한다.

**진행 상황 (0.12.0 당시)** · 로드맵 02를 완료했다. 새 진동수 계산은 단일 IR 그래프를 표시하고, 동위원소 비교 창·CLI는 겹친 스펙트럼과 모드별 수치를 저장한다. 이전 헤시안의 진동수·ZPE 비교도 확인했다. 당시 03–05는 대기 중이었다.

**구현 및 실행 기록** · `engines/pyscf_engine.py`·`frequencies.py`·`infrared.py`와 기존 분자 그림·동위원소 비교 경로를 확장했다. `scripts/validate_infrared.py`의 실계산은 `outputs/ir_validation/20261123_095118_066934_ohtz1m40`에 있다. 여섯 피크의 독립 세기 차이는 최대 0.00130 km/mol, 변위 절반에서 0.000169 km/mol이었다. H₂ 신축은 1.09×10⁻²³ km/mol이었다. WSL 회귀 35건과 숨김 Windows GUI 3건, 실제 CLI의 H₂O/D₂O·HDO·구버전 경로가 통과했다. [사용법·그림·전체 수치](docs/INFRARED_KO.md) · [변경 이력](CHANGELOG.md)

**다음 AI를 위한 인계 (0.12.0 당시)** · [로드맵 03](ROADMAP.md)에서 `gas_mixture.py`의 79/21 몰분율 N₂/O₂ 출력과 `engines/gas_thermo.py`의 종별 NIST 기준 열용량을 같은 1 bar·298.15–2,000 K에서 대조하도록 인계했다.

## QS-20261123-ROADMAP-ISOTOPE · 개발 로드맵과 동위원소 진동 비교

**작업 식별 정보** · ChatGPT Codex · COMPLETED · 2026-11-23 · 0.11.1 · 시작 HEAD `fffb892`

**작업 목적** · 사용자가 요청한 남은 개발 목록을 실행 순서와 완료 기준까지 기록하고 첫 항목을 구현한다. 첫 완료 기준은 저장된 분자 진동수 결과 하나로 H₂O/D₂O 진동수와 영점에너지를 CLI·GUI에서 비교하고, 해석해와 독립 구현에 대조하는 것이다.

**작업 계획**

1. 범용 자체 DFT와 항공우주 물리라는 연구 목표에 맞춰 `ROADMAP.md`에 완료된 기반, 순차 작업, 주제에 따라 선택할 확장을 정리한다. 요구 목록과 현황에서 연결한다.
2. 기존 `engines/frequencies.py`의 강체 운동 투영과 ZPE 계산을 재사용한다. 동위원소 질량·원자 순번 입력, 원본 보존, 별도 결과 저장, `main.py` 및 `molecular_studies_app.py`의 진입점을 추가한다.
3. H₂/D₂의 해석적 질량비, 기존 물 계산의 두 질량 조합과 독립 PySCF 모드, GUI·CLI 경로를 검증한다. 새 결과와 다음 항목을 변경 이력에 남긴다.

**진행 상황** · 로드맵에 순차 과제 5개와 연구 주제별 확장 조건을 작성했다. 01 동위원소 비교의 CLI·GUI·수치 대조를 마쳤다. 02–05와 아래 장기 작업은 남아 있다.

**구현 및 실행 기록** · `molecular_isotopes.py`에서 정규모드와 ZPE를 계산하고, `main.py`·`molecular_studies_app.py`에서 실행한다. 0.8.4 자체 PBE-D4 물 결과의 원본과 배포 예제 모두 D₂O 1176.41/2662.39/2774.07 cm⁻¹을 재현했다. HDO를 포함한 독립 PySCF 대조 및 기울기·저장 등 WSL 회귀 29건, 숨김 Windows GUI 3건을 통과했다. `outputs/isotope_validation/final`·`hdo_final`에 원시 출력이 있다. [수치와 재현 명령](docs/VALIDATION_KO.md#물의-동위원소-진동수-비교) · [변경 이력](CHANGELOG.md)

최종 스테이징 검사에서 Windows로 생성한 예제 JSON의 CRLF 줄끝이 표시되어 LF로 통일했다. 수정 후 두 CLI 치환 결과를 다시 확인했다.

GitHub 자동 push를 수반하는 커밋은 이번 파일의 외부 전송에 대한 명시적 승인을 확인할 수 없다는 자동 승인 검토 결과로 차단되었다. 이 작업의 이력은 커밋 한 번에만 `quantumstarter.autopush=false`를 적용해 로컬에 보존한다.

**다음 AI를 위한 인계 (0.11.1 당시)** · 로드맵 02에서 데카르트 단위 길이의 모드와 IR 계산용 질량 정규화 모드를 구분하고, ±변위 구조의 쌍극자를 수집하도록 인계했다. 당시 배포한 물 헤시안 예제에는 쌍극자 미분이 없다. 이후 0.12.0에서 IR 전용 예제를 추가했다. 아래 별도 장기 HANDOFF_READY 인계가 있어 이 파일을 보존한다.

## QS-20261122-MIX-CONV · 혼합기체 음속과 분자 계산 조건 비교

- 담당 ChatGPT Codex · COMPLETED (2026-11-22, 0.11.0). 시작 HEAD `0091865`. 고정 몰분율 기체의 비열·비열비·음속과 분자 기저·격자 비교를 Windows GUI와 CLI에 구현했다.
- Cantera 3.2.0에서 NASA7 입력 5조성×38온도, 자체 PBE에서 물·H₂ 12조건과 N₂·O₂ 8조건을 검증했다. 2,000 K의 79/21 혼합기체에서 γ=1.4 가정은 +32.12 m/s, 네 DFT 조건 간 최대 차이는 0.190 m/s였다. WSL 회귀 47건·숨김 GUI 2건이 통과했다.
- 원시 계산은 `outputs/mix_conv_validation`에 있다. Cantera 기준 계산은 같은 계수 파일을 UTF-8로 읽고 YAML 문자열로 전달한다. 기저 비교는 상태·모드 순번 기준이며, 분기 개수가 바뀌면 위치를 직접 비교하지 않는다. 이 작업은 완료했고 아래 별도 장기 인계는 보존한다. [재현 명령과 수치 자료](docs/MIXTURE_CONVERGENCE_KO.md) · [변경 이력](CHANGELOG.md)

## QS-20261122-GIT-REMOTE · GitHub 원격과 커밋 후 자동 push

- 담당 Claude · COMPLETED (2026-11-22). 시작 HEAD `65d19ba`, 작업 트리 깨끗함. 사용자 요청: "현재까지 만들어진 것 GIT에 올리고, 편집 끝나면 GIT에 자동 반영". 사용자 선택: Private 저장소, "커밋되면 자동 push"(Claude 응답마다 자동 커밋하는 방식은 고르지 않음).
- 수행: Git Credential Manager로 GitHub 계정 로그인(브라우저 승인), API로 Private 저장소 `mjyoun09/QuantumStarter` 생성, `origin` 연결 후 `main` 커밋 28개 push. `.githooks/post-commit`(자동 push)·`.gitattributes`(훅 LF 고정) 추가, `git config core.hooksPath .githooks`.
- 검증: WSL에서 훅을 직접 실행해 Windows `git.exe` 경유 push 성공(로그 `ok main 65d19ba`). 이 기록의 커밋이 Windows Git Bash 경로의 첫 자동 push 시험이다.
- 다음 AI 주의: 커밋하면 곧바로 GitHub에 올라가므로 비밀값·대용량 산출물을 커밋하지 않는다. 올린 커밋을 amend/rebase하지 않는다(강제 push 안 함). 운영 방법은 PROJECT_STATE의 Git 항목. [변경 이력](CHANGELOG.md)

이전 별도 작업 **QS-20261121-MOLECULAR-SUITE / ChatGPT Codex / COMPLETED (2026-11-22)** — 아래 절.

## QS-20261121-MOLECULAR-SUITE · 분자 계산 3기능 업데이트

- 담당 ChatGPT Codex · COMPLETED (2026-11-22). 시작 HEAD `24518b2`. 자체 단일항 TDA UV–Vis, N₂/O₂/CO₂의 1 bar 열물성·NASA7, 지정 결합의 실수 RKS→UKS 스캔을 분자 GUI와 CLI에 연결했다.
- 최종 원시 계산은 `outputs/molecular_suite/gas_final`과 `outputs/molecular_suite/spectra_scan_final`이다. NIST 비교 108행·NASA7 독립 복원, PySCF와 피크 24개 대조, H₂와 물의 분기, Windows 숨김 GUI와 CLI 실계산, 체크포인트 재시작, WSL 회귀 69건을 확인했다. 첫 반복 TDA 참조의 미수렴은 PySCF 전체 A 행렬 대각화로 해결했다. 재현 시 OMP·OpenBLAS 스레드는 각 1개로 설정한다.
- 이 작업의 추가 수행 항목은 없다. 아래 장기 작업은 HANDOFF_READY로 보존한다. [실행법·검증 자료](docs/MOLECULAR_SUITE_KO.md) · [변경 이력](CHANGELOG.md)

이전 별도 작업 **QS-20261121-H2-DISS / ChatGPT Codex / COMPLETED (2026-11-21)**. H₂ 19개 거리의 자체 RKS·UKS와 PySCF를 대조하고 분기를 1.6178125–1.618125 Å로 좁혔다. 최종 원시 결과는 `outputs/h2_dissociation/20261121_final`, DIIS가 1.70 Å에서 대칭 해로 돌아간 첫 시도는 `20261121`에 보존했다. 이 자료는 0.10.0의 일반 분자 실수 RKS→UKS 스캔을 검증할 때 사용했다. [검증과 재현 명령](docs/VALIDATION_KO.md) · [변경 이력](CHANGELOG.md)

이전 별도 작업 **QS-20261121-N2-THERMO / ChatGPT Codex / COMPLETED (2026-11-21)**. 질소 1 bar, 298.15–2,000 K의 PBE/RRHO 36점과 NIST 비교, 좌표 변위 절반 민감도 검사, 열역학 회귀 4건을 마쳤다. 원시 계산은 `outputs/n2_thermo/20261121`와 `20261121_halfstep_verified`에 있다. 앞선 `20261121_halfstep`는 중단된 계산이다. [표와 재현 명령](docs/VALIDATION_KO.md) · [변경 이력](CHANGELOG.md)

이전 별도 작업 **QS-20261121-ATOM-MIX / ChatGPT Codex / COMPLETED (2026-11-21, 0.9.8)**. 원자 SCF의 Anderson 혼합을 구현했다. H–Kr PBE 72쌍 및 대표 4종×3회 측정이 모두 통과했다. PySCF·Libxc 포함 WSL 회귀 62건 통과. 장기 작업은 아래 HANDOFF_READY 상태로 보존한다. [검증 수치와 재현 명령](docs/VALIDATION_KO.md) · [변경 이력](CHANGELOG.md)

이전 별도 작업 **QS-20261117-ATOM-GGA / ChatGPT Codex / COMPLETED (2026-11-18, 0.9.7)**. H–Kr 비분극·분극 PBE 72건, 닫힌 껍질 7종과 붕소의 독립 에너지 대조, 기존 LDA/LSDA 72건 원본 및 NIST 회귀를 완료했다. [방법과 재현 명령](docs/ATOMIC_PBE_KO.md) · [검증 수치](docs/VALIDATION_KO.md) · [변경 이력](CHANGELOG.md)

## 작업 식별 정보

- 작업 ID: QS-20261029-NATIVE-DFT
- 작업명: 첨부 범용 DFT 요구사항의 누락 구현, 자체 SCF 기반부터 단계별 검증
- 담당 AI: Claude — 2026-11-01 사용자 "codex가 12단계 완료하고 커밋 했다는데 해볼래"(직전 대화에서 Claude가 제안한 IRC·열역학 독립 검증 진행 승인)에 따라 단계 13 인수. 시작 HEAD 2dd6c8f, 작업 트리 깨끗함. 이전: ChatGPT Codex(단계 10·12).
- 상태: HANDOFF_READY — 단계 13 COMPLETED(2026-11-01, CHANGELOG 0.9.4). 0.10.0에서 실수 RKS→UKS 결합 스캔과 세 기체의 열물성·NASA7을 구현했다. 장기 요구에는 quasi-RRHO·용매, 복소 궤도 안정성, 항공우주 세부 연구 주제 선정이 남아 인계서를 보존한다. **GUI 시험은 QUANTUMSTARTER_TK_HIDDEN=1로만.**
- 별도 작업 QS-20261112-HE-LDA(중성 He 방사 격자 LDA-SCF, 0.9.5, Claude)는 2026-11-12 완료 — CHANGELOG 참조. 작업 메모 TASK_HE_LDA_SCF.md는 사용자 지시에 따라 완료 후 삭제했다.
- 별도 작업 QS-20261117-ATOM-LSDA(원자 LDA/LSDA H–Kr, 로그 방사 격자, 0.9.6, Claude)는 2026-11-17 완료 — CHANGELOG 참조. 작업 메모 TASK_ATOM_LSDA.md는 완료 후 삭제했다. 이 방향의 다음 후보는 아래 "방사 격자 원자 DFT 다음 단계" 절.

### 원자 GGA(PBE) 최초 계획 (2026-11-17 제안, 2026-11-18 완료)

사용자가 0.9.6 완료 후 기록을 요청한 최초 계획을 보존한다. 완료 근거는 이 문서 첫 줄 및 CHANGELOG 0.9.7에 있다.

- 출발점: HEAD `1bec0db`(0.9.6). 엔진 radial_atom.py, 검증 scripts/validate_atoms.py, 시험 test_radial_atom.py. he-lda(atomic_dft.py)는 0.9.5 기록 보존을 위해 건드리지 않는다.
- 목표: 같은 로그 격자 엔진에 PBE(Libxc `GGA_X_PBE` 101 + `GGA_C_PBE` 130, PySCF 문자열 `PBE` = (101, 130) — 2026-11-17 parse_xc로 확인)를 LDA/LSDA와 나란히 추가. 스핀 비분극·분극 모두. 런타임 Libxc 없이 닫힌 식으로 구현하는 원칙 유지.
- 구현 요점
  1. 교환: F_x(s) = 1 + κ − κ/(1 + μs²/κ), κ = 0.804, μ ≈ 0.2195, s = |∇n|/(2k_F n). 스핀 분극은 E_x[n↑, n↓] = ½(E_x[2n↑] + E_x[2n↓]).
  2. 상관: PW92 LDA 상관 + H(rs, ζ, t). GGA_C_PBE가 내부에서 쓰는 PW92 변형(`LDA_C_PW` 12 또는 `LDA_C_PW_MOD` 13)과 β·γ 상수 자릿수는 구현 전에 Libxc 값 대조로 확정한다(기억으로 정하지 않는다).
  3. 구대칭 밀도에서 ∇n = n′(r)r̂, σ = n′². 퍼텐셜 v = ∂e/∂n − (1/r²) d/dr[r²·2(∂e/∂σ)n′](스핀은 σ↑↑, σ↑↓, σ↓↓). 로그 격자에서 d/dr = (1/J)d/dx. 이산 E_xc의 정확한 미분이 되도록 차분을 고르지 않으면 이중 계산 식과 직접 식이 어긋나므로, 두 식의 차이를 반드시 기록한다(현재 LDA/LSDA는 ≤1e-14).
  4. 주의: 원점(첨점, n′ ≠ 0)과 바깥 꼬리(n → 0에서 s 발산, F_x는 1 + κ로 포화)의 수치 잡음, Libxc 밀도 문턱값(1e-15)과의 꼬리 차이(0.9.6 발견). LDA/LSDA 경로의 결과는 비트 단위로 바뀌면 안 된다(72건 회귀).
  5. CLI: `--xc`는 분자 옵션과 겹치므로 원자용은 별도 이름(예: `--functional lda|pbe`)으로 하고 설정 JSON에 저장한다.
- 검증 계획(NIST 표에는 GGA가 없다)
  1. 점별: e_xc, vrho, vsigma를 Libxc(PBE, spin 0/1)와 상대 ~1e-12, 퍼텐셜은 ∂E/∂n 수치 미분과 대조.
  2. 닫힌 껍질 He·Be·Ne·Mg·Ar·(Zn·Kr): PySCF RKS `PBE`, 0.9.6과 같은 짝수 템퍼드 s(+p, +d는 Zn·Kr) 기저 수열, 기준 |ΔE| ≤ 1e-6 Ha(Richardson).
  3. 열린 껍질: 구면 평균을 맞추려면 PySCF `scf.addons.frac_occ`(존재 확인됨) 같은 분수 점유가 필요하다. 쓸 수 있는지·같은 밀도를 주는지 먼저 확인하고, 안 되면 열린 껍질은 Libxc 점별 대조와 격자 수렴만으로 검증하고 한계를 적는다.
  4. 문헌 PBE 원자 에너지표는 출처를 원문에서 확인한 경우에만 쓴다.
  5. 격자 수렴 차수 2, 전자 수, 미수렴 처리, 기존 LDA/LSDA 72건 NIST 대조 재통과.
- 당시 후보 중 Anderson 혼합은 0.9.8에서 완료했다. 남은 후보는 비구면·다중항 구분, 스칼라 상대론(ZORA/Koelling–Harmon), 원자 결과를 이용한 노름 보존 의사퍼텐셜 생성이다. 순서는 미정이다.

### 단계 14 당시 재개 후보

아래는 0.9.4 당시의 계획이다. 이후 NASA7과 실수 스핀 해리 검증을 완료했다. 현재 작업 순서는 [개발 로드맵](ROADMAP.md)에서 관리한다.

- 사용자 응용 방향: DFT × 항공우주 물리. 구체 주제는 사용자·지도교수 결정 대기. 열물성 방향이 확정되면: 온도 구간별 Cp/H/S 표 출력, NASA 7계수 다항식 적합(200–1000/1000–6000 K 관례), JANAF/NASA Glenn/ATcT 등 기준 데이터 대조(출처·라이선스 확인 후). 반응 속도 방향이면: TS·IRC(0.9.4) 위에 TST k(T)(분배함수 비, 터널링 보정은 별도 결정).
- 코드 측 남은 항목: quasi-RRHO·용매/표준상태, 삼중항(RKS→UKS)·복소 SCF 안정성. IRC는 단계 13 완료이므로 반복 구현하지 않는다.
- 재개 시: PROJECT_STATE → 이 문서 → `git status --short`, `git log -5 --oneline`. 최근 검사 로그 `outputs/irc_validation/wsl_final.log`(83개 통과). 이미 통과한 회귀는 이유 없이 반복하지 않는다.

### 단계 13: 분자 IRC + 열역학 독립 교차검증 (Claude, 완료 — CHANGELOG 0.9.4)

- 목적: TS에서 양 방향 최소 에너지 경로(IRC)를 따라 반응물·생성물 연결을 확인하는 `calculation='irc'`. 항공우주 방향 후보(TST 속도)의 전제인 "TS가 원하는 두 최소점을 잇는가"를 확인하는 도구. 함께 RRHO 열역학을 PySCF 외 독립 구현(ASE IdealGasThermo)과 교차검증.
- 결정(2026-11-01, Claude): IRC 구동기는 **geomeTRIC 1.1.1 내장 IRC**. 근거: 설치 확인 결과 Sella·jax 미설치(Sella는 jax 추가 의존), geomeTRIC IRC는 기존 TS 탐색과 같은 내부좌표·엔진 어댑터를 쓰고 Gonzalez–Schlegel(1990) 질량가중 2단계 스텝, 'both' 방향, 끝 근처에서 최소화 전환을 제공(소스 optimize.py IRC_step/evaluate_IRC_step/reset_irc 확인). geomeTRIC 질량은 자체 PeriodicTable(표준 원자량)이므로 경로 정의 질량을 결과에 기록. 거부된 스텝은 progress에서 제거됨(optimize.py:981) → 경로는 geomeTRIC progress를 사용. 단일 분자 TRIC는 run_optimizer가 DLC로 전환(optimize.py:1291). ASE Calculator 래퍼는 이번 범위에서 불필요(열역학 비교는 ASE 함수에 진동수만 전달).
- 완료 기준: (1) NH₃ 반전 TS에서 IRC 양방향 수렴, 두 끝점 에너지 = 별도 최적화 NH₃(허용오차는 최적화 기준 수준), 두 끝점은 거울상, 경로 에너지가 TS에서 양쪽으로 단조 감소. (2) HCN⇌HNC: TS 탐색→허수 진동 1개→IRC→끝점이 독립 최적화 HCN·HNC와 구조(RMSD)·에너지 일치 — 검증 스크립트·outputs/ 보존. (3) 최소점 입력은 명시적 오류, 제약과 함께 요청은 입력 오류. (4) ASE IdealGasThermo와 H·S·G 일치(비선형·선형·원자). (5) 관련 회귀 + 숨김 GUI 통과. (6) 문서·CHANGELOG·버전 0.9.4.
- 계획: (a) engines/optimizer.py의 geomeTRIC 엔진 어댑터를 공유 함수로 분리(동작 불변) → (b) `irc_molecule`: run_optimizer(irc=True, 'both', hessian='first'), 경로·질량가중 호 길이 s·상대에너지·양 끝점·장벽·단조성·끝점 RMSD(고유/거울 포함), TS 재계산을 실행 폴더 최상위에 → (c) molecular.py CALCULATIONS·제약 규칙, engines/__init__.py, main.py, molecular_app.py(콤보 설명·지표) → (d) test_irc.py, test_frequencies에 ASE 비교 → (e) scripts/validate_irc.py HCN/HNC → (f) 회귀·문서.
- 변경 예정: engines/optimizer.py, molecular.py, engines/__init__.py, main.py, molecular_app.py, test_irc.py(신규), test_frequencies.py, scripts/validate_irc.py(신규), README_KO·docs(VALIDATION/REQUIREMENTS_STATUS/NATIVE_DFT)·PROJECT_STATE·CHANGELOG, 버전 문자열.
- 예상 문제: 선형 끝점(HCN/HNC)에서 내부좌표 특이점, IRC 끝 근처 진동(geomeTRIC이 최소화로 전환), 미수렴 시 run_optimizer 예외로 progress 소실 → 자체 기록으로 진단.
- 진행(완료):
  - (a) optimizer.py: `_recording_engine`/`_start`/`_geometric_molecule`로 분리(동작 불변; 오류 문구만 "… 계산을 중단합니다"). (b) 신규 engines/irc.py `irc_molecule` + `superposition_distance`/`rmsd`(가중 Kabsch, 거울상 허용 선택). (c) molecular.py CALCULATIONS에 'irc'(원자 ≥2, 제약 불가는 기존 규칙), engines/__init__.py, main.py, molecular_app.py(설명·`_irc_metric` 지표). (d) test_irc.py(입력·Kabsch·NH₃ 양방향·최소점 거부·라이브러리 부재), test_frequencies `test_ideal_gas_thermo_matches_ase`, test_molecular_gui `test_real_native_irc_run`(wait_finished timeout 인자 추가). (e) scripts/validate_irc.py(HCN⇌HNC).
  - 발견·수정: geomeTRIC 1.1.1 run_optimizer는 IRC일 때 `M.molecules`(조각 수)를 좌표계 생성 전에 읽는데, customengine 분자는 토폴로지가 없어 AttributeError → irc.py에서 `molecule.build_topology()` 선호출(파일 입력 경로와 같은 거리 기반 결합).
  - 실측(탐색, scratchpad): NH₃ native PBE/STO-3G/grid3, TS 탐색 29.3초 후 IRC 41.3초. 양방향 수렴, 경로 33점, 끝점 ΔE 0.0011 kJ/mol, 거울상 RMSD 7.6e-5 Å(고유 회전 RMSD 0.443 Å), 오르막 0, s=±1.625 bohr·amu^½, 끝점 장벽 54.99 kJ/mol, TS 허수 −880.6 cm⁻¹ 1개, 기울기 계산 63회(Hessian 24). 거부된 스텝은 경로에서 제외 확인.
  - Windows: test_irc 입력 2개 + FrequencyInputTests 4개 통과. ASE 비교 첫 실행에서 물 298 K의 G만 1.93e-7 eV 차로 상대 2e-6 기준 초과 → G=H−TS는 큰 항의 차이라 허용오차 기준을 |H|+T·S로 변경(허용 상대값 2e-6 유지). 측정: S 상대차 ≈4.0e-7 = ASE kB(CODATA2014)와 3.39e-7 차 + 소량, H ≤2.9e-7. 물 S(298 K) 188.56 J/mol/K.
  - HCN⇌HNC 첫 실행(outputs/irc_validation/hcn_hnc_pbe_svp_20261101/): 연결은 맞았으나 HNC(선형) 끝에서 geomeTRIC IRC가 23스텝 후 ~1e-6 Å 스텝으로 77스텝 정체 → 방향 한도 100 도달 후 geomeTRIC은 조용히 반대 방향으로 넘어가 "Converged" 보고. 우리 끝점 검사가 HNC 끝(최대 기울기 9.9e-4 > tight 4.5e-5)을 미수렴으로 잡음. 조치: 기준 미달 끝점 일반 최소화(endpoint_<방향>/), `irc_end`·`relaxation` 기록, 마지막 방향 한도 시 geometric_irc.xyz(평가마다 기록, Å 10자리)를 기록된 계산과 맞춰 경로 복구(state['gradients'] 추가), 성공 판정 = 두 끝점 기준 충족. test_irc에 방향당 12단계 한도 시험 추가.
  - HCN⇌HNC 재실행(…_relax/): IRC 수렴 181.6초, forward 103(한도)·backward 34회, HNC 끝 최소화 4단계(이동 RMSD 4.2e-4 Å, ΔE −6.4e-7 Ha) 후 독립 HNC와 RMSD 1.29e-4 Å·ΔE 2.0e-8 Ha, HCN 끝 RMSD 5.9e-5 Å·ΔE −5.4e-8 Ha, 끝점 진동수 허수 없음, TS 허수 −1077.13 cm⁻¹ 1개. 장벽 194.08/132.39, 이성질화 61.68 kJ/mol(PBE/def2-SVP; 문헌 대조 안 함).
  - CLI(Windows→WSL) NH₃: outputs/irc_validation/cli/20261101_144237_018561_debc906r, exit 0, 0.9.4, 탐색 결과와 같은 값. CLI 콘솔에 OPT_JSON 줄이 출력되는 것은 기존 optimize와 같은 동작(변경하지 않음).
  - 검사: WSL 최종 83개 통과(166.8초), test_irc 6개(58.3초), 리팩터 직후 회귀 56개(74.4초), Windows 숨김 실계산 GUI 3개(72.8초), Windows 빠른 시험 46개(19 통과/27 skip). git diff --check·py_compile·변경 파일 LF 바이트 검사 통과. 문서: README_KO, VALIDATION 0.9.4, REQUIREMENTS_STATUS, PROJECT_STATE, CHANGELOG.

### 단계 12 완료 기록 (ChatGPT Codex, 2026-10-30 구현·검증 / 2026-11-01 정리)

- 목적/완료 기준: 기존 FFTDF bulk 힘·응력을 사용한 원자/셀 동시 최적화, 명시적 에너지(E 또는 smearing F)·힘·응력 수렴 판정 및 반복 한도/실패 구분. 고정 정수 FFT mesh·scaled k점을 유지하고 모든 단계/최종 구조를 보존한다. FFT mesh별 E/F·힘·응력 차이를 비교하는 보조 기능과 GUI/CLI 설정을 제공한다.
- 계획: (1) ASE 셀 필터/최적화 API와 응력 단위 검토, (2) PeriodicConfig·새 periodic optimizer·dispatch·CLI, (3) mesh 수렴 도구, (4) periodic_app 설정·결과 및 숨김 GUI 시험, (5) 독립 에너지 최소점/미분·실계산·실패/취소·회귀 검증, (6) 문서/CHANGELOG/Git 보존.
- 주요 변경: periodic.py, engines/periodic.py, engines/periodic_optimizer.py(신규), engines/__init__.py, main.py, periodic_app.py, periodic_convergence.py(신규), dft_workflow.py/dft_worker.py, energy_analysis.py/periodic_plotting.py, 셀/mesh/GUI 검증 파일 및 scripts/, 관련 Markdown.
- 완료: 셀 완화·FFT 검사·CLI/GUI·취소/최종 구조 비교 구현. Si 등방(7단계, F)/전단 full(10단계, E) 실제 수렴과 인근 체적 에너지 상승, 초기·최종 FFT21³/25³/29³ 검사 통과. 숨김 GUI6개 통과. 상세 방법 문서 작성.
- 검증 완료: WSL 고유51개 + 숨김 GUI6개 통과. core_regression27개/65.286초, remaining_tests8개/162.337초, 공유 실행17개/29.430초(자동 mesh1개 중복), gui_final6개/26.300초. `outputs/cell_validation/`에 파일 로그와 반환 코드 보존(공유 실행17개는 도구 출력). 첫 밴드 회귀 묶음은 요약 없이 exit1로 종료했으나 로그 runner로 분리 재검사해 모두 exit0. 원인은 확정 못 했으며 허용오차를 완화하지 않았다.
- 최종 Windows CLI 확인: `main.py --config examples/silicon_cell.json --no-show --output outputs/cell_validation/final_cli`, exit0. `20261030_202319_940271_g2izg_ts/`: version0.9.3, 7단계 수렴, 힘1.44886e-4 eV/Å·평균응력 잔차7.48862e-4 GPa. `settings.json=cell_optimize`, `final_settings.json=gradient`, `final/` 재검사, 루트 최적화 판정·궤적·전자구조·periodic.png 보존 확인. 그림은 저장 파일로 확인했고 GUI 창은 띄우지 않음.
- 진행 중: 없음. Python 구문 컴파일·git diff --check 통과. 결과 데이터는 Git 제외이므로 별도로 보존한다. 최종 코드/문서 커밋 제목은 `Add bulk cell relaxation and FFT mesh convergence workflows`.
- 미수행: 전체 장기 요구 중 IRC·quasi-RRHO·용매·삼중항/복소 안정성 및 항공우주 구체 연구 주제는 이후 항목으로 보존한다.
- 다음 재개: 이 절을 우선한다. 새 최적화의 SCF 미수렴을 성공 처리하지 않으며, coarse FFT mesh의 수치 수렴을 물성 정확도로 설명하지 않는다.

### 단계 13 재개 후보 (Codex 작성 — 단계 13에서 IRC로 진행·완료됨, 역사 기록)

- 사용자 응용 방향(PROJECT_STATE): HEA가 아니라 **DFT × 항공우주 물리**. 구체 주제는 사용자·지도교수 결정 대기. 제안된 1순위: 고온 기체·추진제의 온도별 열물성(Cp/H/S, NASA 7계수 다항식)과 TST 반응 속도, JANAF/NASA 기준 독립 검증 — 현재 엔진(경원소 분자, TS/NEB, 진동수·RRHO)과 가장 잘 맞음.
- 코드 측 남은 항목: IRC, quasi-RRHO·용매/표준상태, 삼중항·복소 SCF 안정성. 주기계 셀 완화·FFT 검사·GUI 노출은 단계12 완료이므로 반복 구현하지 않는다.
- 우선 후보 IRC: `engines/optimizer.py:optimize_molecule`의 geomeTRIC Engine 어댑터, `engines/frequencies.py:frequency_molecule`/`harmonic_analysis`, `test_constraints.py`의 NH₃ TS 검증을 읽는다. 설치 geomeTRIC `params.py`에 `irc`/`irc_direction` 인자가 존재함만 확인했으며 실제 동작·질량 가중 경로·양쪽 종료 판정은 아직 검증하지 않았다. 구현 전 ACTIVE_TASK에 별도 단계 목표·완료 기준을 적는다. 허수 모드1개 확인된 TS에서 양 방향 추적·최소점 최적화·연결 확인을 검증 후보로 삼고, 실행 완료를 반응 연결의 증명으로 처리하지 않는다.
- 재개 시: PROJECT_STATE → 이 문서 → `git status --short`, `git log -5 --oneline` → 해당 모듈. `Get-Content outputs/cell_validation/core_regression/tests.json -Encoding UTF8`로 최근 검사 확인. 이미 통과한 회귀는 이유 없이 반복하지 않는다.
- 현재 제한: bulk/FFTDF/외부압력0만 셀 완화, 원자 제약·최적화 재시작·포논 안정성·자동 k점/기저 수렴은 없음. 실험 조건/실패 기록은 docs/VALIDATION_KO.md 0.9.3에 있다.

### 단계 11: bulk FFTDF 힘·응력 API 통합 (Claude, 완료 — CHANGELOG 0.9.2)

- 계획: (1) scripts/validate_periodic_extensions.py `fftdf()` 확인, (2) 일반 k점(고정 scaled k)·UKS·smearing에서 PySCF FFTDF 기울기/응력 지원 여부와 FD 일치 실측, 지원 안 되면 설정에서 명시 거부, (3) PeriodicConfig 명시적 `integral_scheme`(gdf 기본/fftdf)·계산 종류, FFT mesh 기록, (4) engines/periodic.py 분기·힘(eV/Å)·응력(GPa) 단위·부호·출처, (5) CLI·energy_analysis 호환(적분 방식 혼합 금지), (6) 테스트·문서. GDF+Becke 에너지와 FFTDF+uniform 힘/응력을 섞지 않는다. 셀 최적화는 그 다음.
- 실측(Si 변위, GTH-SZV, PBE, FFT mesh 25³, UniformGrids, scratchpad/probe_fftdf.py): UKS(spin 0) 기울기 vs FD 7.4e-9 Ha/bohr, 응력 xx vs 변형 FD 9e-12. smearing(Fermi 0.1 eV) 기울기·응력은 **자유에너지 e_free의 미분**과 일치(1e-8, 1e-11), 내부에너지 e_tot 미분과는 1.8e-4 차 → 결과에 '미분 대상 = F'를 명시해야 함. 2×2×2 k 메쉬: 기울기 vs FD 4.6e-9, 응력 xx vs FD 3.3e-11 → 일반 k점·UKS·smearing 모두 지원 확인.
- 구현: PeriodicConfig `integral_scheme`(gdf/fftdf), `fft_mesh`, `calculation`(single_point/gradient) + 규칙(fftdf는 bulk만, gradient는 fftdf만, fft_mesh는 fftdf만). engines/periodic.py: FFTDF면 density_fit 생략·UniformGrids(같은 mesh)·small_rho_cutoff 0, `periodic_gradient`(힘 eV/Å, 응력 Ha/bohr³·GPa, 압력 −tr/3, 미분 대상 E 또는 F 명시, mesh·부피 기록), 결과 integral_scheme/fft_mesh/xc_grid. energy_analysis 주기계 호환에 integral_scheme·fft_mesh 추가. CLI --integral-scheme, --fft-mesh, 주기계 --calculation gradient 허용.
- test_periodic_gradient 4개 통과(61.7초; 첫 실행 입력 시험 1개는 내 시험 사례가 calculation=gradient를 유지해 다른 규칙이 먼저 걸린 시험 오류 → 사례 수정): 제품 API로 기울기 vs FD(1e-7), 응력 xx vs 변형 FD(1e-9), GPa·eV/Å 변환, 압력, smearing은 F 미분, GDF/FFTDF 반응 분석 혼합 거부.

### 단계 10 완료 및 최신 인계 (ChatGPT Codex, 2026-10-30)

- 목적/완료 기준: 기록된 Cr·Co 정체 재현 → 원인 수정 → 정수 점유, 기존 잔차/밀도/에너지 기준을 유지한 수렴 및 실수 내부 안정성 확인. Newton 응답을 독립 유한차분·PySCF 기준과 비교하고 기존 RKS/UKS 회귀를 통과한다. 다음으로 bulk FFTDF 힘·응력의 실제 지원 범위와 slab 밴드를 확인한다. IRC/quasi-RRHO/용매 등 장기 항목은 별도 범위를 명시하며 미완료를 숨기지 않는다.
- 계획: (1) 재현/계측을 outputs/에 보존, (2) native_newton/native_stability 및 필요 시 자체 XC 선형응답 모듈 수정, (3) test_newton·응답 검증 추가와 관련 회귀, (4) 주기계 지원 검토·slab 계산, (5) 상태/검증/요구사항/CHANGELOG 갱신.
- 변경 예정: engines/native_newton.py, engines/native_stability.py, 필요 시 engines/native_response.py, test_newton.py 및 응답 검증 파일, scripts/ 검증 스크립트, 관련 Markdown.
- 완료: Cr·Co 미수렴 재현과 수정, 해석적 응답 검증·성능 계측, 기존 물·Fe·OH 회귀, slab 실계산, bulk FFTDF 도입 가능성 검토, 관련 문서 및 완료 이력 갱신. 과거 완료 기록은 유지.
- 진행 중: 없음. 코드·계산 결과 저장 완료.
- 재현: WSL 기본 샌드박스 접근은 E_ACCESSDENIED; 승인된 `wsl --exec` 경로로 실행. `scripts/validate_newton.py --case Cr --max-cycle 50 --label baseline_20261030`: 기존 코드 144.8초, 잔차 3.9276e-5, 미수렴(격자3/PBE/STO-3G/스핀6/1스레드). outputs/newton_validation/에 이력 보존.
- 구현 중: native_response.py 자체 J/K·XC 선형응답(Libxc Cartesian f_xc만 외부), native_stability 해석적 H·v/기존 FD 비교 경로, native_newton 궤도 간격 전처리·CG 정밀도·반올림 크기 수락 검사. 0/1차원 안정성 고유문제 처리.
- 검증: test_native_response 8조건(RKS/UKS × LDA/PBE/PBE0/B3LYP)의 응답·Hessian FD, PySCF 응답, 대칭성/선형성, 캐시 없는 스트리밍 통과(5.45초). 첫 시험의 core 궤도 극단 스핀 밀도에서 FD 오차 ~1e-6 발생; 비정상점에 작은 일반 궤도 회전을 주어 노드/스핀 포화의 비매끄러움을 피한 뒤 허용오차 완화 없이 통과. 임의 가산 밀도 대신 물리적 궤도 회전으로 변위.
- Cr/Co 실계산 성공: outputs/newton_validation/analytic_pcg_20261030/{Cr,Co}/report.json. Cr 73회/117.37초, E=-1032.8088678851316, 잔차2.10e-13, 최저곡률9.51e-7. Co(스핀3) 52회/64.68초, E=-1367.1350701366355, 잔차1.03e-12, 최저곡률2.36e-4. 기존 물리/수렴 기준·정수 점유 유지. 서로 다른 사례 동시 실행 시간이라 엄격 성능 벤치마크로 일반화하지 않는다.
- OH 회귀 첫 실행: 수렴/에너지/최종 안정성은 통과하나 명시적 안장점 탈출 분기 기록을 요구하는 시험 단언 실패. 전처리 Newton 자체가 음의 곡률을 따라 탈출 가능하므로 구현 경로 단언을 최종 안정성 확인으로 교체. 실패를 숨기지 않으며 전체 재검증 예정.
- 첫 관련 회귀: WSL 54개 통과(174.029초), Cr·Co 최종 밀도에서 독립 PySCF Newton도 수렴하며 에너지 차1e-9 Ha 이내. He/STO-3G의 회전 차원0 경고를 발견하여 Steihaug의 0 잔차 조기 종료 수정. 0/1차원·음의 곡률·최종 정준화 후 수렴 기준까지 추가 검사하는 최종 회귀 진행 중.
- 속도 계측: 물 PBE/STO-3G/grid3/1스레드, 같은 수렴 궤도·방향·각7회 H·v 중앙값 해석0.032245초 vs FD0.120587초(3.74배), 최대 차3.04e-8. 초기 캐시 구축 포함 첫 응답0.07681초. 전체 SCF 속도 보장으로 일반화하지 않는다. outputs/newton_validation/benchmark_20261030/water/report.json.
- 주기계 검토: scripts/validate_periodic_extensions.py. graphene(GTH-SZV, 2×2×1, 16 Å 진공, PBE, GMKG16점) 수렴, kz=0, 메쉬 고유값 차3.80e-9 eV, DOS 적분 오차3.31e-7 states. 결과 outputs/periodic_extension_validation/stage10_20261030/slab/.
- FFTDF bulk 검토: 변위 Si, gamma/PBE/GTH-SZV, uniform FFT 25³/스미어링 없음. 기울기 vs SCF 에너지 FD 차7.39e-9 Ha/bohr, xx 응력 vs 고정 분율 좌표 변형 FD 차9.27e-12 Ha/bohr³. xy는 거의0이라 비영(非零) yz 성분도 추가 확인 예정. 기존 GDF/API는 바꾸지 않았으며 FFTDF가 별도 bulk 경로로 가능하다는 검토 결과다(일반 k점/UKS/smearing/메쉬 수렴·입력/출력/GUI 통합은 후속).
- 추가 확인: Fe 자체 Newton 53회 수렴(E=-1249.5569908940256, 잔차7.24e-13, 최저곡률4.31e-4). Fe H·v도 FD 대비3.74배(각7회 중앙값), 최대차1.62e-8. FFTDF 비영 yz 응력=0.0003806094765 Ha/bohr³, SCF 에너지 FD 차1.25e-12. 각각 benchmark_20261030/Fe, stage10_shear_20261030/fftdf에 보존.
- 최종 검증 중 주의: `-W error`로 실행한23개에서 native/응답/안정성10개와 주기계5개 통과, ASE `atoms.py`의 NumPy2.5 shape 할당 DeprecationWarning 때문에8개 오류. 제품 수치 오류가 아니며 환경 변경 없이 해당8개만 `-W error::RuntimeWarning`으로 재검증 중. 실패 원문 final_focus_wsl.log 보존.
- Windows 첫21개는 로그 OK이나 PowerShell stderr 리다이렉트가 프로세스 성공과 다른 종료코드1을 반환. `outputs/newton_validation/run_windows_checks.py`가 subprocess 실제 종료코드를 수집하도록 재검증 중(환경 숨김 모드 유지).
- 최종 결과: 위 ASE 경고로 중단된8개를 수치 RuntimeWarning=오류로 재검증하여 모두 통과(24.802초). WSL 고유69개 확인(첫54개 + 새 수치 경계2개 + 주기계13개; 중복 실행은 합계에 제외). 최종 native/Newton/응답/안정성10개는 모든 경고=오류에서도 통과. Windows21개 재검증77.997초, 실제 returncode=0, skip 없음(`windows_verified.log`, `windows_exit_status.json`). Python 컴파일·git diff --check·변경 파일 LF 바이트 검사 통과. 외부 ASE deprecation은 알려진 환경 경고로 남기며 패키지/환경은 수정하지 않았다.
- 완료 이력: CHANGELOG 0.9.1 및 docs/VALIDATION_KO.md에 조건·실측·실패/조치·재현 명령 보존. 출력 데이터는 .gitignore 대상이므로 outputs/는 현 PC에서 별도로 보존해야 한다.
- 미수행(다음 작업): bulk FFTDF 힘·응력의 제품 API 통합, IRC, quasi-RRHO·용매, 삼중항·복소 안정성. HEA 정밀 기준값은 사용자 현재 응용 목표가 아니므로 자동 우선순위로 삼지 않는다.
- 정확한 재개 지점: bulk를 확장한다면 `scripts/validate_periodic_extensions.py`의 `fftdf()`를 먼저 읽는다. 일반 k점(고정 scaled k), UKS·smearing 자유에너지 미분 및 mesh 수렴을 검증한 뒤 `PeriodicConfig`의 명시적 FFTDF/균일 격자 설정, `engines/periodic.py` 분기·힘/응력 단위/출처, CLI/GUI·energy_analysis 호환성 검사를 추가한다. 현재 GDF+Becke 에너지를 FFTDF+uniform 힘/응력과 섞지 않는다. 셀 최적화는 그 다음이다.
- 우선 명령: `git status --short`; `Get-Content docs/VALIDATION_KO.md -Encoding UTF8 -TotalCount 45`; `wsl --exec .venv-dft/bin/python scripts/validate_periodic_extensions.py fftdf --label <새이름>`. 이미 통과한 Cr·Co/전체 회귀는 변경 이유 없이 반복하지 않는다.
- 아래 단계 2~9 및 말미의 첫 단계 인계는 역사 기록이며 현재 미구현 목록으로 해석하지 않는다. 이후 작업을 시작하면 이 최신 절의 담당자·상태·계획부터 갱신한다.

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

### 단계 10 인수 당시 재개 목록 (수행 결과는 위 최신 인계 참조)

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

## 단계 12 중간 기록 (2026-10-30 당시 기록, 최종 결과는 위 완료 절)

- `PeriodicConfig.cell_optimize`, ASE BFGS/FrechetCellFilter 기반 원자·셀 완화, 고정 FFT 정수 격자/상대 k점, E 또는 smearing F 목적함수, 힘·응력·에너지 개별 종료 기준 구현.
- `periodic_convergence.py`: 고정 구조 mesh 3개 이상, 마지막 연속 두 구간의 E/F·힘·응력 차이 판정, JSON/CSV 및 개별 계산 보존. CLI 및 GUI 연결 작성 중.
- 첫 Γ-only Si / FFT25³ / smearing 0 검증은 step 3 SCF 미수렴으로 명시적 중단. `outputs/cell_validation/first_si` 보존. 미수렴 힘을 최적화에 전달하지 않음.
- smearing 0.1 Γ-only 시도는 지나친 체적 팽창 및 느린 SCF를 보임. 이를 물리적 Si 평형 결과로 해석하지 않음. 예제 kmesh 2³, conv_tol 1e-9로 변경하여 별도 검증 예정.
- PowerShell `python -c` 문자열 인자 따옴표 소실/Store alias 실행 실패는 검증 harness 문제. 파일 스크립트 및 `.venv/Scripts/python.exe`, WSL `.venv-dft/bin/python` 사용.
- 구현 초안이며 아래 후속 단계는 검증 전: 실제 셀 수렴/인근 체적 E/F 확인, mesh 실제 실행, 실패·취소·최대단계 회귀, 숨김 GUI, 문서와 버전 갱신.

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

0.8 첫 단계 종료 당시의 다음 구현 지점은 아래와 같다. 최신 순서는 [개발 로드맵](ROADMAP.md)에 있다.

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
