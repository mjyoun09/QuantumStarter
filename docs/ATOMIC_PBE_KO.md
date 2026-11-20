# 로그 방사 격자 원자 PBE · 0.9.7

작업 QS-20261117-ATOM-GGA, 2026-11-17~18, ChatGPT Codex. `atomic_pbe.py`와 `radial_atom.py`는 기존 NumPy/SciPy 환경에서 원자 PBE를 계산한다. 검증 스크립트에서 PySCF 및 Libxc와 대조한다.

## 실행과 범위

```powershell
.\.venv\Scripts\python.exe main.py --model atom --element Ar --functional pbe --no-show
.\.venv\Scripts\python.exe main.py --model atom --element N --functional pbe --spin-polarized --no-show
.\.venv\Scripts\python.exe main.py --model atom --config <결과폴더>\settings.json --no-show
```

`--functional lda|pbe`는 원자 전용이며 기본값 `lda`는 기존 Slater+VWN5와 같다. 분자 `--xc`와 혼용하지 않는다. 선택값은 settings/results JSON에 저장되고 콘솔·그림·XC 식별자에 반영된다. 수렴/미수렴/오류 종료 코드는 0/2/1이다.

H–Kr, 비상대론, 고정 점핵, 구면 평균 밀도, 껍질별 고정 점유만 다룬다. 열린 껍질에서 각 m 성분을 균등하게 채우고 스핀 분극 시 다수 스핀 우선이다. 기본 배치는 기존 NIST 배치를 유지하며 PBE로 배치 간 바닥상태를 탐색하지 않는다. 다중항·비구면·상대론·음이온·혼성 함수는 이 원자 엔진의 범위 밖이다.

### Anderson 밀도 혼합 (0.9.8)

```powershell
.\.venv\Scripts\python.exe main.py --model atom --element Fe --functional pbe --spin-polarized --mixing-method anderson --no-show
```

`--mixing-method`를 생략하면 기존 설정 파일도 선형 혼합으로 실행된다. `--mixing`은 두 방법의 감쇠 계수 β(기본 0.3)다. 초기 밀도, 고정 점유와 격자, SCF 종료 조건은 공통이다. 매번 실제 출력 밀도에 대해 |ΔE| < 10⁻¹⁰ Ha와 혼합 전 스핀별 적분 절대 잔차의 합 < 10⁻⁷ e를 함께 검사한다(기본 설정).

이력 5개의 감쇠 Type-II Anderson 방법을 쓴다. 입력 방사 밀도를 x, 출력과 입력의 차이를 R, 선형 혼합 결과를 g = x + βR라 두면 새 입력은 g − Δg·γ다. Δ는 연속 반복 사이의 차이를 열로 모은 행렬이다. γ는 ∑σ∫|Rσ − (ΔR·γ)σ|² dr를 최소화하며 한 계수 벡터를 두 스핀에 공통으로 쓴다. 최소제곱에 로그 격자의 적분 가중치 J·h를 반영한다. [Walker와 Ni의 고정점 Anderson 가속 방법](https://epubs.siam.org/doi/10.1137/10078356X)

첫 갱신은 선형이다. 최소제곱 풀이 실패, 열 종속(상대 특잇값 절단 10⁻¹⁰), 계수 절댓값 합 > 20, 유한하지 않은 후보 또는 어느 스핀이든 적분 음전하 크기가 그 스핀 전자 수의 10⁻⁶을 넘으면 차분 이력을 비우고 해당 갱신에 선형 혼합을 쓴다. 허용 범위 안의 음수는 0으로 자르고 각 점유 스핀을 초기 전자 수에 맞춰 정규화한다. Anderson 실행의 `scf_history.csv`와 `results.json`에서 `input_mixing_step`을 보면 각 반복의 입력을 만든 방법을 확인할 수 있다. [H–Kr 비교 및 재현 명령](VALIDATION_KO.md)

## 함수와 미분

교환은 PBE의 `F_x(s)=1+κ−κ/(1+μs²/κ)`, `s²=σ/[4(3π²)^(2/3)n^(8/3)]`이다. 분극 교환은 스핀 스케일링 `E_x[↑,↓]=(E_x[2↑]+E_x[2↓])/2`를 쓴다. 상관은 PW92 수정형과 PBE H의 합이다. 근거: [PBE 원 논문](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865), [Libxc 7의 PBE 정의](https://gitlab.com/libxc/libxc/-/raw/7.0.0/maple/gga_exc/gga_c_pbe.mpl).

상수는 κ=0.804, β=0.06672455060314922, μ=βπ²/3, γ=(1−ln2)/π². PW92는 일반 ID12가 아니라 ID13의 A=(0.0310907,0.01554535,0.0168869)와 f″(0)=1.7099209341613656을 사용한다. `GGA_X_PBE` 101 + `GGA_C_PBE` 130이며 LDA 옵션의 VWN5와 다르다. [PW92 수정형 상수 및 구분](https://gitlab.com/libxc/libxc/-/raw/7.0.0/src/lda_c_pw.c), [PBE 공통 상수](https://gitlab.com/libxc/libxc/-/raw/7.0.0/maple/util.mpl).

지역 에너지 `e=n ε`를 정의하고 2개(비분극 n,σ) 또는 5개(분극 n↑,n↓,σ↑↑,σ↑↓,σ↓↓) 변수에 자체 1차 chain rule을 적용한다. 유한차분이나 외부 자동미분을 런타임에 쓰지 않는다. 상관 꼬리의 ε+H 상쇄는 다음 동치식으로 피한다:

`ε_c(PBE)=γφ³ log1p[expm1(y)/(1+u+u²)]`, `y=ε_c(PW92)/(γφ³)`, `u=A t²`.

총밀도 ≤1e−30 bohr⁻³에서는 0을 반환한다. 스핀 비율 `2nσ/n`의 분수 거듭제곱에 1e−30 바닥을 적용하고 같은 에너지를 미분한다. 정확히 빈 스핀에서 소수 스핀 PBE 미분은 특이하므로 그 퍼텐셜은 수치 규약이다. 점유된 스핀과 에너지 극한은 유지한다. Libxc의 1e−15 밀도 절단과 다르며 진공/완전 분극 끝점에서 동일한 수치를 약속하지 않는다.

## 방사 격자의 변분 일관성

기존 로그 격자 `r=a(exp(x)−1)`, `J=r+a`, `a=1e−4/Z`를 쓴다. `D n=J⁻¹ d n/dx`는 내부 중앙차분, 양 끝 2차 한쪽 차분이다. 원점에서 n′=0을 강제하지 않는다. 밀도 적분 가중치는 `W=4πr²Jh`.

`E_xc=Σ_i W_i e(n_i,(Dn)_i)`에 대해 `v_xc=e_n+W⁻¹DᵀW e_grad`를 사용한다. 구면 발산의 이산형이며 끝점 항까지 포함한 정확한 전치이다. 분극 flux는 `2 e_σ↑↑ n′↑+e_σ↑↓ n′↓`와 그 반대 스핀이다. 다른 차분으로 −div를 다시 근사하면 같은 이산 에너지의 미분이 되지 않는다.

시험은 전체 격자의 임의 방향·매끄러운 방향 에너지 유한차분, 끝점 포함 전치 항등식, SCF 직접 에너지와 이중 계산식, 격자 2차 수렴을 각각 확인한다.

## 검증 범위와 재현

- WSL Python 3.12.3, NumPy 2.5.4, PySCF 2.14.0, Libxc 7.0.0. Windows Python 3.13.5. 환경·패키지 변경 없음.
- 일반 스핀 비율 0.01..0.99, 밀도 1e−6..1e3의 180점/스핀 모드에서 개별 x/c 및 합의 ε·vrho·vsigma 대조. 상대 기준3e−12; 상쇄되는 합 vsigma는 |x|+|c| 크기로 정규화한다. σ=0에서 PW92_MOD 일치도 검사한다.
- 극단 비율은 Libxc의 `1±ζ` 계산에서 작은 스핀의 유효 숫자를 잃는다. 별도 60자리 Decimal의 원래 ε+H 식, 상대 변위1e−18 중앙차분으로 6점의 x/c 모든 미분을 대조해3e−12 통과했다. 일반 대조 허용오차를 넓혀 처리하지 않았다.
- H–Kr PBE 72건(비분극/분극), 기본 h=.001 모두 수렴. He·Ar·Kr 비분극, H·N·Fe 분극에서 h=.002/.001/.0005 관측 차수2. 전자 수·마디 수·직접/이중계산 에너지 검사 포함. 원자별 원시 출력과 이력을 보존한다.
- 기존 LDA/LSDA 72건은 HEAD26381a7 원본 소스와 현재 소스를 **동일 Windows 프로세스**에서 실행했다. 에너지·성분·고윳값·전체 SCF 이력·밀도·퍼텐셜·파동함수가 모두 비트 단위로 같다. 과거 WSL JSON과의 마지막 비트 차이는 같은 환경의 원본 비교로 구분했다. NIST 기본 에너지 오차도 기존 범위다.
- 열린 껍질 B(2p¹)는 PySCF UKS의 `frac_occ_`로 독립 검증했다. 두 기저에서 다수 스핀 2p 점유는 (⅓,⅓,⅓), 소수 스핀은 0이었다. SCF 매 반복과 최종 밀도를 반지름 25개·방향 5개에서 비교해 최대 각도별 밀도 차이 7.21×10⁻¹¹ e/bohr³를 확인했다. 방사 Richardson 에너지와의 차이는 1.83×10⁻⁷ Ha였다. 다른 열린 껍질의 비구면·다중항 상태에는 이 결과를 일반화하지 않는다.
- Zn·Kr의 독립 기준은 별도 구대칭 가우시안 RKS 계산으로 검증했다. He·Ar을 함께 계산해 PySCF RKS와 비교 가능한 조건을 확보했다. Zn·Kr의 약형 에너지 및 기저 수렴 절차는 `scripts/validate_spherical_gaussians.py`에 있다.

```powershell
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_atomic_pbe.py pbe --label <새이름>
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_atomic_pbe.py reference --label <같은이름>
.\.venv\Scripts\python.exe scripts/validate_atomic_pbe.py compare --label <같은이름>
.\.venv\Scripts\python.exe scripts/validate_atomic_pbe.py lda-source --label <새이름>
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_atomic_pbe.py fractional --label <새이름>
wsl --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .venv-dft/bin/python scripts/validate_spherical_gaussians.py --radial-report outputs/atomic_pbe_validation/<같은이름>/pbe/report.json --out outputs/atomic_pbe_validation/<새폴더>
wsl --exec .venv-dft/bin/python -m unittest test_atomic_pbe test_radial_atom test_atomic_dft test_solver test_workflow -v
```

`pbe --resume`는 같은 폴더의 완료 원자부터 건너뛴다. 계산 조건/수식이 바뀌었다면 새 label을 사용해야 한다. `reference --heavy`는 선택적 Zn/Kr 기준을 추가한다. `lda`는 과거 `validate_atoms` 보고서 대조, `lda-source`는 Git 원본과의 엄격한 비트 비교다. `outputs/`는 Git 제외이므로 원시 결과는 현 PC에 별도 보존한다.

## 발견과 조치

1. Cr 분극 꼬리에서 매우 작은 소수 스핀 기울기의 제곱은 언더플로하고 교차 곱은 남을 수 있다. 최초 Gram 입력 검사에서 이를 비물리 입력으로 잘못 거부했다. 입력 허용 검사에만 float tiny를 써 수정했으며 실제 에너지/미분 수식을 바꾸지 않았다. 전용 시험 추가 후 Cr부터 재개했다.
2. PySCF 2.14는 기본 겹침 고유값 절단1e−6으로 조밀한 even-tempered 기저 일부를 제거했다. 독립 기준에서는1e−10으로 명시하고 대각화와 DIIS가 같은 직교 공간을 쓰게 했다. 과도한 지수 상한1e9는 매우 높은 가상 궤도 때문에 낮은 고유상태의 정밀도를 떨어뜨렸다. 원자에 맞는 상한과 기저 수열로 수렴을 확인한다. 400→1000→2000 방사 적분점 변화는 Be에서 약1e−12 이하라 첫 수 µHa 차이의 주원인이 아니었다.
3. 원래 기준 탐색의 미수렴 결과는 기준으로 쓰지 않는다. 중단한 일반 수축 기저/과도한 조밀 기저 실험과 로그도 삭제하지 않았다. 자체 원자 SCF 수렴 기준은 변경하지 않았다.

원자별 외삽 에너지 오차와 기준 계산의 수렴 조건은 [검증 기록 0.9.7](VALIDATION_KO.md)에 있다.
