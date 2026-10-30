# 자체 분자 DFT · 0.9.1 구현과 검증

첨부 프롬프트에서 가장 먼저 요구한 자체 SCF를 추가했다. 기존 PySCF 백엔드를 유지하며, 같은 화면의
**계산 방법 → 계산 엔진**에서 `native` 또는 `pyscf`를 선택한다. 기본값은 기존 결과와의 연속성을 위해
`pyscf`다. **자체 DFT · 물 (LDA 검증)** 예제는 native로 설정된다.

## 직접 구현한 부분과 라이브러리 경계

| 부분 | 구현 |
| --- | --- |
| 점유수, 밀도행렬, 실공간 밀도/밀도 기울기 | engines/native_scf.py |
| J/K의 ERI 수축, LDA/GGA XC 행렬 조립, 하이브리드 정확교환 계수 | 자체 코드 |
| 겹침행렬 직교화, 일반화 KS 고유문제 구성, SCF, Pulay DIIS, 밀도 혼합 | 자체 코드; 행렬 고유값/선형풀이에는 NumPy/SciPy 사용 |
| 총에너지 성분 합산, 에너지/밀도/commutator 수렴 판정 | 자체 코드 |
| Newton·실수 내부 안정성·해석적 궤도 Hessian 응답 | native_newton.py / native_stability.py / native_response.py |
| 기저 데이터, 겹침/운동/핵인력/전자반발 적분, AO 값과 도함수 | PySCF/libcint |
| 원자 중심 방사·각도 격자 및 분할 가중치 | PySCF |
| XC 에너지 밀도와 1·2차 밀도 미분 | Libxc (PySCF의 libxc 인터페이스; 2차 응답은 Cartesian 밀도/gradient 미분) |
| CUBE 내보내기와 별도 밀도 적분 교차검사 | 기존 PySCF 내보내기 도구 |
| pyscf 백엔드의 SCF 전체 | 외부 PySCF; 자체 구현 성과로 계산하지 않음 |

native는 PySCF의 RKS/UKS `.kernel()`, `get_veff()`, DIIS, 밀도/XC 행렬 조립을 사용하지 않는다.
외부 kernel 호출을 오류로 막은 상태에서도 자체 실계산이 통과하는 테스트가 있다.
공통 분자 구성/파일 출력은 기존 engines/pyscf_engine.py에 있으며 실제 SCF 드라이버가 명시적으로 분기된다.

## 밀도와 에너지 규약

AO 겹침행렬을 S, MO 계수를 C라 하면 CᵀSC=I다.
RKS는 점유수0 또는2, D=C n Cᵀ, N=Tr(DS). UKS는 스핀별 점유수0 또는1, D=Dα+Dβ다.

J[pq]=Σrs (pq|rs)D[rs], Kσ[pq]=Σrs (pr|qs)Dσ[rs].
RKS F=H+J+Vxc−aK/2, UKS Fσ=H+J+Vxc,σ−aKσ이며 a는 Libxc의 global hybrid 계수다.

E = Tr(DT) + Tr(DVnuc) + Tr(DJ)/2 + Exc,semilocal + Ex,exact + Enuc.
정확교환 항은 RKS −a Tr(DK)/4, UKS −a ΣσTr(DσKσ)/2다.
오비탈 에너지 합을 총에너지로 쓰지 않는다. 저장 성분의 합이 총에너지와1e-7 Ha 이상 다르면 실패한다.

LDA는 `lda,vwn`: Slater 교환과 VWN 상관이다. 결과의 `provenance.xc_functional_ids_and_weights`,
`exact_exchange_fraction`, `environment.libxc`에 실제 Libxc 함수 ID/계수/버전을 기록한다.
PBE에는 실제 밀도 기울기를, PBE0/B3LYP에는 실제 K 행렬을 사용한다. 선택 D4 분산 보정(0.8.2)은 dftd4 라이브러리가 계산하는
밀도 무관 기하 보정이며 SCF가 끝난 뒤 별도 성분 `dispersion`으로 더한다. 자체 SCF 알고리즘에는 관여하지 않는다.

SCF는 에너지 변화, 밀도행렬 Frobenius 변화, 직교 기저에서 FDS−SDF의 Frobenius 잔차를 **모두** 만족해야 한다.
최종 정수 점유·canonical 궤도도 같은 기준으로 다시 확인하며 느슨한 추가 수렴 기준을 쓰지 않는다.
미수렴에서는 마지막 궤도 시도의 밀도/에너지를 진단값으로 저장한다. 중간 혼합밀도의 이력과 최종 진단값은 다를 수 있다.

## 현재 범위와 제한

- native: H~Kr, 지원 기저의 비상대론 전체전자·점핵, RKS/UKS, LDA/VWN·PBE·PBE0·B3LYP.
  원소별 기저(`basis_overrides`, 0.8.1)도 같은 AO 기저 구성으로 사용한다. 혼합 기저 물에서 외부 엔진과 3.1e-11 Ha 일치(VALIDATION 0.8.1).
- 직접 dense ERI를 보유하므로 최대160 AO이며 메모리 추정치가 `memory_mb`를 넘으면 실행을 거부한다. 모든 H~Kr 화합물의 수렴/정확도를 검증한 것은 아니다.
- 기존 pyscf: H~Og의 Dyall 기저, ECP, X2C 등 기존 지원 유지. native에서 미지원인 옵션은 오류로 명시한다.
- 자체 Newton·구조 최적화·진동·NEB는 구현됐다. 자체 ECP/X2C 및 주기계 자체 SCF는 미구현이다. 분산 에너지·기울기는 D4(외부 dftd4)로 지원한다.
- **해석적 원자 기울기(0.8.3, engines/native_gradient.py)**: dE/dR_A = ΣD·dh/dR_A + J' + K' + XC' − ΣW·dS/dR_A + dV_nn/dR_A.
  W = ½DFD(RKS) 또는 Σ_s D_sF_sD_s(UKS), 2전자 미분 `int2e_ip1`은 한 bra 셸씩(메모리 3·n_shell·n³),
  XC는 고정 격자에서 ρ와 ∇ρ의 핵 좌표 미분(AO 1·2차 도함수)으로 계산한다. 격자 가중치 미분은 없다.
  고정 격자 유한차분과 1.2e-9~3.3e-9 Ha/bohr, 같은 격자의 PySCF 기울기와 1e-12 이하로 일치(VALIDATION 0.8.3).
- **DIIS 수정(0.8.4)**: Pulay 계수 상한 30이 수렴 근처의 정상적인 큰 계수를 막아 OH에서 300회 중 289회 이력을 버리고 정체(ΔD ~2e-5)했다. B 행렬 고유분해 의사역행렬(수치적 영공간만 제거)로 바꿔 OH가 ~25회에 수렴한다. 다른 계의 반복 수·에너지는 불변.
- **수렴 ≠ 최소**: 그 OH 해는 PySCF 안정성 분석에서 내부 불안정(더 낮은 UKS 해 존재, ~8e-6 Ha)으로 판정됐다. 근축퇴 열린 껍질은 여러 SCF 해가 있을 수 있다. 이후 추가한 자체 실수 내부 안정성 분석으로 검출할 수 있다.
- 구조 최적화에서 단계당 자체 SCF 반복 수가 PySCF보다 많았다(물 PBE-D4 예: 8~14회 vs 5~7회). 0.8.6에서 개선(아래).
- **초기 추정·혼합(0.8.6)**: `auto` 초기 추정 = 원자 밀도 중첩(PySCF `init_guess_by_minao`, 외부 보조). 투영 밀도는 전자 수가 정확하지 않으므로 첫 Fock 생성에만 쓰고 그 점유 궤도에서 시작한다. 밀도 혼합은 교환자 노름 > `native_damping_residual`(기본 1e-3)일 때만. 물 14→9, CH4 14→7, NH2 14→10, O2 10→8회, Fe 원자(STO-3G, 2S=4) 수렴. 과거 남았던 Cr·Co 사례는 아래 0.9.1 Newton으로 해결됐다.
- **SCF 안정성 분석(engines/native_stability.py, `native_stability=check|follow`)**: 궤도 회전 C·exp(K)에 대한 에너지 기울기 g_ai = 2n·F_ai(n=2 RKS, 1 UKS)를 직접 구성한다. 0.9.1부터 Hessian-벡터 곱은 해석적 응답, 최저 고유값은 Lanczos(eigsh; 0/1차원은 별도 처리)다. 기존 중심차분은 `hessian_vector_fd`로 검증용 보존한다. follow는 불안정 모드 방향 선 탐색 후 SCF를 다시 수행(최대3회). **실수 내부 회전만 검사**하며 RKS→UKS(삼중항)·복소 불안정은 검사하지 않는다.
- **Newton 개선(0.9.1)**: 고정 격자에서 δF=J[δD]+f_xc[δD]−aK[δD](RKS 교환 계수1/2)를 직접 조립하고 궤도 회전 항을 더해 H·v를 계산한다. Libxc `eval_xc_eff(deriv=2)`는 밀도/Cartesian 기울기에 대한 함수 도함수만 제공한다([공식 인터페이스](https://pyscf.org/pyscf_api_docs/pyscf.dft.html)). AO·f_xc를 메모리 예산 내 캐시하고 부족하면 블록 단위로 재계산한다. 궤도 에너지 간격 전처리(최솟값1e-3 Ha), 잔차에 따른 CG 정밀도, 에너지 반올림 수준에서 기울기 감소를 요구하는 수락 검사로 기존 정체를 해소했다. 최종 정준화 뒤 잔차·에너지·밀도 기준을 다시 확인한다. Cr·Co(PBE/STO-3G/grid3/스핀6·3) 수렴·양의 실수 내부 곡률 확인. 가상의 분수 점유나 느슨한 수렴 기준은 사용하지 않았다.
- 자체 프로그램이라는 것은 모든 수치 라이브러리를 새로 썼다는 뜻이 아니다. 적분/XC 라이브러리는 두 비교 엔진이 공유하므로 비교는 주로 자체 SCF/행렬/에너지 구현을 검증한다.

## 중단, 체크포인트, 재시작

SCF 반복마다 `checkpoint.npz`를 임시 파일에 쓴 뒤 원자적으로 교체한다. pickle을 사용하지 않는다.
중단 버튼은 `cancel.request`를 만들며 native의 적분 블록/SCF 또는 PySCF SCF 콜백에서 반영된다.
큰 외부 적분·SCF 한 단계 중에는 요청 반영까지 기다릴 수 있다. 프로세스를 강제로 죽이지 않는다.
첫 반복 전 중단이면 체크포인트가 없을 수 있다. `cancelled`는 성공·미수렴·실패와 별개이며 CLI 종료코드는3이다.

실행 기록에서 체크포인트 재시작을 선택하거나 CLI `--restart`를 쓴다. 기본 버튼은 과거 설정을 불러와 실행한다.
물리 조건·원소 순서·기저·전하/스핀·격자 지문과 밀도 크기·유한성·대칭성·스핀별 전자 수·점유수 범위를 검사한다.
반복 한도/수렴 기준/자원/DIIS 전략 변경은 허용한다. 같은 물리 조건이면 native↔pyscf 간 밀도를 전달할 수 있고
출처 엔진을 기록한다. DIIS 이력/정확한 실행 시점을 복원하는 것이 아니라 저장 밀도에서 새 SCF를 시작한다.
개발 중 schema1 체크포인트와 schema2는 구분하며 호환되지 않는 버전은 오류로 거부한다.

## 사용성 추가

- 계산 전 3차원 구조 미리보기와1부터 시작하는 원자 번호로 거리/각도 측정.
- 빠른 확인/정밀 격자 프리셋: 바뀌는 grid_level, conv_tol, 반복 한도, 그림 크기, 자체 잔차 기준을 화면에 표시.
- SCF 실시간 탭: 에너지, 밀도 변화, 잔차. native commutator와 PySCF orbital gradient는 정의가 다르다고 명시.
- 결과 JSON의 `energy_components_hartree`에서 성분 확인. 실행 기록의 엔진 출처, 다중 선택 비교, 흡착에너지 분석.

## 실행

```powershell
.\.venv\Scripts\python.exe main.py --config examples\native_water.json --engine native --no-show
.\.venv\Scripts\python.exe main.py --config examples\native_water.json --restart <기존실행폴더> --no-show
wsl --exec .venv-dft/bin/python -m unittest test_native test_energy_analysis test_molecular test_periodic -v
wsl --exec .venv-dft/bin/python scripts/validate_native.py
$env:QUANTUMSTARTER_GUI_DFT_TEST='1'
$env:QUANTUMSTARTER_TK_HIDDEN='1'
.\.venv\Scripts\python.exe -m unittest test_molecular_gui -v
```

## 실제 검증 결과

최종 회귀: WSL 자체/분석/기존 분자·주기계39개, Windows 실계산 GUI10개(skip 없음), 기존 수소 workflow5개 통과.
원본 로그와 테스트 범위는 [검증 기록](VALIDATION_KO.md)에 보존했다.

영구 비교: `outputs/native_validation/20261029_195245_027252/comparison.json`.
공통 STO-3G / grid_level2 / conv_tol1e-10 / 2스레드. H·He·H₂·H₂O는 LDA/VWN, OH는 PBE다.
에너지는 Ha이며 근사·기저의 연구 정확도를 뜻하지 않는다.

| 계 | 자체 에너지 | 독립 PySCF 에너지와 차이 | 자체 밀도에서 PySCF 재시작 후 차이 |
| --- | ---: | ---: | ---: |
| H (2S=1) | -0.435670232979561 | -1.11e-16 | -1.11e-16 |
| He | -2.771886044438531 | -4.44e-16 | -4.44e-16 |
| H₂ | -1.121206106145701 | 4.44e-16 | 4.44e-16 |
| H₂O | -74.732103087545160 | -7.11e-14 | -4.26e-14 |
| OH (2S=1) | -74.572249584501690 | -2.40551e-6 | -2.84e-14 |

OH 독립 실행의 차이를 숨기거나 허용오차를 느슨하게 해서 같은 상태라고 처리하지 않았다.
같은 밀도에서 평가하면 에너지 차이4.26e-14 Ha, Fock 최대차3.55e-15 Ha이며,
외부 SCF를 그 밀도에서 시작해도 같은 에너지로 수렴한다. 퇴화 열린껍질 궤도/각도 격자 및 SCF 해 선택의
영향으로 해석할 수 있지만 모든 스핀·궤도 해를 탐색해 원인을 확정한 것은 아니다.

추가 테스트: RKS/UKS×4 XC의 Fock가 에너지의 유한차분 미분과2e-7 이내 일치,
9조건의 기준 비교, CᵀSC·D 재구성·Tr(DS)·스핀·성분 합, 초기 추정(core/seed1729 perturbation),
격자 변경, 강제 미수렴, 중단·재시작과 불일치 설정 거부. 초기 추정 비교는 바닥상태 완전 탐색이 아니다.

공식 인터페이스 참고: [PySCF DFT](https://pyscf.org/user/dft.html),
[Libxc 인터페이스](https://pyscf.org/_modules/pyscf/dft/libxc.html),
[AO 값/도함수 API](https://pyscf.org/_modules/pyscf/dft/numint.html).
