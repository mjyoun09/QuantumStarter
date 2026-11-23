# 첨부 범용 DFT 프롬프트 구현 추적

기준: 2026-10-29 사용자 첨부 프롬프트. 상태는 실제 코드와 실행 증거로 갱신한다.
프롬프트의 우선 구현 대상은 자체 분자 DFT다. 앞으로 만들 기능의 순서·선행 조건·완료 기준은 [개발 로드맵](../ROADMAP.md)에 모았다.

| 요구 | 현재 상태 | 근거/후속 작업 |
| --- | --- | --- |
| 수소 검증기 보존, 입력/GUI/엔진/저장 분리 | 구현 | solver.py, molecular.py, dft_workflow.py, engines/, 별도 GUI |
| 원소·좌표·단위·전하·스핀·격자·PBC | 구현 | 분자 외부 엔진 118원소, 주기계는 기저별 지원. 분자 구조 최적화에서 원자 고정·거리 제약을 지원한다. |
| 자체 분자 DFT: 밀도/Fock/SCF/DIIS/에너지 | **구현·실계산 검증** | native_scf.py; RKS/UKS, LDA/GGA/global hybrid. 적분·XC 함수는 외부 라이브러리 |
| 실공간 구대칭 원자 Kohn–Sham | **H–Kr LDA/LSDA 및 비분극·분극 PBE 구현·검증** | `--model atom --functional lda\|pbe`. 원자 PBE 자체 수식, 이산 에너지 미분, 72건 수렴, 닫힌 껍질 7종 독립 대조 ≤1e−6 Ha. LDA/LSDA 72건 원본과 동일. 비상대론·점핵·고정 배치·구면 평균. [상세](ATOMIC_PBE_KO.md) |
| RKS/UKS, LDA/GGA/hybrid, ECP/X2C, 밀도/MO | 자체/외부에서 구현, 범위 다름 | 자체 ECP/X2C는 미지원. 외부 PySCF의 기존 지원 보존 |
| 에너지 성분, 초기 추정 선택·비교 | 구현·검증 | 성분 합·전자 수·잔차·core/perturbed 검증 |
| 중단·체크포인트·재시작 | 분자 구현·검증, 주기계 중단 구현 | 분자 물리 설정/전자 수 검사·재시작. 주기계는 SCF/셀 단계 경계 취소와 이력 보존, 최적화 재시작은 미지원 |
| 실시간 SCF, 응답 GUI, 기록/내보내기 | 구현·검증 | 잔차 정의/엔진 구분, 수렴 그래프, 다중 기록 비교 |
| 구조 미리보기·거리/각도, 프리셋 | 구현·검증 | geometry_tools.py, 실제 변경 수치 표시 |
| 원소별 기저 선택 | **구현·검증 (분자)** | `basis_overrides`, GUI/CLI/JSON. 원소별 ECP·X2C 판정, 반응 분석의 원소별 대조. 주기계는 미지원 |
| 분산 보정 | **분자 D4 에너지·힘 구현·검증** | 선택 의존성 dftd4 4.3.1. PBE/PBE0/B3LYP, 110원소, 절단 없는 분자 합. 주기계 D4는 후속 과제다. |
| DFT 원자 힘 (분자) | **구현·검증** | 자체 해석적 RKS/UKS 기울기(LDA/GGA/global hybrid), PySCF 경로, D4 기울기. 격자 응답 항 없음 |
| 구조 최적화 (분자) | **구현·검증** | geomeTRIC(선택 의존성) + 자체/PySCF 해석적 힘. PySCF 자체 geomeTRIC 드라이버와 구조 일치 |
| 고정/거리 제약·스캔 (분자) | **구현·검증** | geomeTRIC `$freeze`/`$set` 제약, 이완 스캔(점마다 제약 최적화). 고정은 geomeTRIC 허용오차(~1e-4 Å) 내 |
| Hessian·진동·ZPE (분자) | **구현·검증** | 해석적 힘 유한차분 Hessian, 정확한 병진·회전 제거, 허수 진동 판정, ZPE·진동 열역학. PySCF 해석적 Hessian과 ≤0.5 cm⁻¹ |
| 동위원소 진동수·ZPE 비교 | **구현·검증 (0.11.1)** | 저장된 헤시안과 구조에 전체 H→D 또는 원자별 질량을 적용한다. H₂/HD/D₂ 해석해, H₂O/HDO/D₂O PySCF 대조. [수치와 재현 명령](VALIDATION_KO.md#물의-동위원소-진동수-비교) |
| 분자 IR 적분 세기·동위원소 스펙트럼 | **구현·검증 (0.12.0)** | 고정 원점의 쌍극자 미분, 질량 정규화, 원본·치환체 겹침 그림과 가변 가우스 선폭. 물 6개 피크, H₂ 비활성, 변위 절반과 독립 방향 미분 검증. [그림과 수치](INFRARED_KO.md) |
| SCF 안정성 분석 (자체) | **구현·검증** | 내부 회전 check/follow. 결합 거리 스캔에서는 실수 RKS→UKS 곡률·분기·하위 해를 추적. 복소 궤도 검사는 미구현. [검증](MOLECULAR_SUITE_KO.md) |
| 닫힌껍질 UV–Vis, 기체 NASA7 | **구현·검증 (0.10.0)** | 단일항 TDA의 에너지·세기와 N₂/O₂/CO₂의 1 bar 열물성·두 구간 계수. [사용법과 대조표](MOLECULAR_SUITE_KO.md) |
| 고정 조성 혼합기체 및 분자 기저·격자 비교 | **구현·검증 (0.11.0)** | NASA7 비열비·음속을 Cantera와 대조. 물·H₂ 및 N₂/O₂ 기체 계수의 기저·격자 변화량과 실행 시간 기록. [수치·재현 명령](MIXTURE_CONVERGENCE_KO.md) |
| 자체 2차(Newton) SCF | **구현·검증** | 신뢰영역 Newton + 안장점 탈출. 0.9.1 해석적 J/K/XC 응답·전처리·수락 검사 개선, FD는 검증용 보존 |
| 전이금속 원자 SCF (자체) | **기록된 Fe·Cr·Co 사례 수렴** | 0.9.1 Cr(스핀6)·Co(스핀3)/PBE/STO-3G/grid3 정수 점유 수렴, 실수 내부 안정성·독립 PySCF Newton 확인. 임의 화합물/스핀 탐색 보장 아님 |
| 병진·회전 열역학, 자유에너지 G | **구현·검증 (분자, RRHO)** | ZPE/H/G 결과 필드, 대칭수 자동/입력, 반응 분석 선택. PySCF thermo와 H·G 1e-10 상대. 0.9.4: 독립 구현 ASE IdealGasThermo와 5종(비선형·선형·원자, 단일·이중·삼중항) × 298/1500/3000 K에서 H ≤3.4e-7·S ≤4.1e-7 상대(ASE의 구 CODATA k_B 차 3.4e-7 수준). quasi-RRHO·용매·표준상태 미구현 |
| 반응/흡착에너지 및 호환성 검사 | 구현·검증 (제한 있음) | energy_analysis.py: 수렴/조건/수지/출처 검사. 주기계 동일 셀·k점만 허용, ZPE/열역학 없음 |
| 전이상태 탐색 (분자) | **구현·검증** | geomeTRIC P-RFO(초기 Hessian은 우리 힘의 유한차분). NH₃ 반전 TS: 허수 진동 정확히 1개 확인 |
| NEB (분자) | **구현·검증** | geomeTRIC CI-NEB + 자체/PySCF 힘, 선형 보간 초기 경로, 양 끝 고정. NH₃ 반전: 대칭 경로, 최고 이미지 = TS 탐색 장벽 ±0.001 kJ/mol |
| IRC (분자) | **구현·검증 (0.9.4)** | `calculation=irc`: geomeTRIC IRC(Gonzalez–Schlegel) + 자체/PySCF 힘, 양방향, 끝점별 기울기 판정·기준 미달 끝점 최소화, 단계 한도 시 경로 복구. NH₃ 반전(거울상 끝점), HCN⇌HNC(끝점 = 독립 최적화 이성질체, 진동수 최소점). 주기계·제약·재시작 없음 |
| 주기계 k점·smearing·2D/3D 정전기 | 외부 PySCF에서 구현 | GDF+Gaussian/GTH, 평면파 아님 |
| 결정·표면·진공/합금 후보 비교 | 부분 구현 | HEA/Pt 등 fcc111 후보 생성 |
| 주기계 힘·응력 (bulk) | **구현·검증 (0.9.2 API, 0.9.3 GUI)** | `integral_scheme=fftdf` + `calculation=gradient`. Γ·일반 k점·UKS·smearing(F 미분)을 FD로 확인. slab·GDF는 미지원(PySCF 2.14). 응력 일부는 내부 적분 유한차분 |
| 주기계 셀 최적화 (bulk) | **구현·검증 (0.9.3)** | `cell_optimize`, ASE BFGS/FrechetCellFilter, full/isotropic, 고정 FFT·상대 k점, E/F·힘·응력 판정. Si 등방 및 전단 변형 완화 검증 |
| FFT mesh 수렴 보조 | **구현·검증 (0.9.3)** | 고정 구조, 마지막 두 연속 구간 E/F·힘·응력 비교, CLI/GUI·개별 실행 보존. k점/기저 수렴은 별도 |
| 표면 전하 보정 | 미구현 | 주기계는 중성 셀만 지원 |
| 밴드·DOS (주기계) | **구현·검증 (bulk 및 slab)** | DOS 항상(가우시안), 밴드는 band_path 지정 시(ASE 경로, get_bands 비자기일관). Si 및 graphene 2D GMKG 경로 실계산 검증 |
| 자체/외부/ML 출처 구분 | 구현 | 결과 SCF 드라이버·적분/XC 라이브러리·체크포인트 출처 구분. 다른 엔진 에너지 혼합 차단 |
| 용매·전체 선형응답 TDDFT·MD | 향후 필요 시 | 닫힌껍질 단일항은 TDA 근사를 지원 |

모든 항목이 끝난 것으로 표시하지 않는다. 실제 검증 결과는 VALIDATION_KO.md 및 별도 자체 엔진 검증 문서에 기록한다.
