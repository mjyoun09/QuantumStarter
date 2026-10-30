# QuantumStarter 0.8 · 자체/외부 분자 DFT와 표면 탐색

**자체 SCF 엔진을 추가했습니다.** `계산 방법 → 계산 엔진 → native` 또는 **자체 DFT · 물 (LDA 검증)** 예제를 선택하세요.
밀도·Fock/XC 행렬·DIIS·SCF·에너지 성분은 직접 구현하고, 기저/적분/격자는 PySCF, XC 함수는 Libxc를 사용합니다.
자체는 H~Kr 전체전자·비상대론·점핵이며 기존 PySCF 118원소 경로와 구분됩니다. [자체 엔진과 실제 검증](docs/NATIVE_DFT_KO.md).
중단/체크포인트/재시작, 실시간 SCF 그래프, 구조 미리보기·거리/각도, 수치 프리셋도 추가했습니다.
**실행 기록**에서 조건별 결과 비교와 [반응·흡착에너지 분석](docs/ENERGY_ANALYSIS_KO.md)을 사용할 수 있습니다.
첨부 프롬프트 전체의 후속 연구 기능이 끝난 것은 아닙니다. [미구현 항목 추적](docs/REQUIREMENTS_STATUS_KO.md)을 확인하세요.

**분자 계산을 H~Og, 118개 원소로 확장했습니다.** Dyall-v2z/v3z, 스칼라 X2C, 유한 핵 모형을 추가했습니다.
상단 **원소별 지원표**에서 사용 가능한 기저와 주기계 지원을 확인하세요. [원소 지원 안내](docs/ELEMENTS_KO.md).
주기계와 ML 모델의 원소 범위는 별도이며, 모든 원소의 고체·표면 계산을 뜻하지 않습니다.

개발 AI는 [AGENTS.md](AGENTS.md) → [PROJECT_STATE.md](PROJECT_STATE.md) → ACTIVE_TASK.md(존재 시)를 먼저 읽습니다.

**고체·표면 계산을 추가했습니다.** 상단 **고체·표면 / GPU 탐색** 버튼에서 격자·k점·2D/3D 경계 조건을 입력합니다.
Ag·Au·Cu·Pd·Pt 표면의 H 후보를 MLIP로 빠르게 탐색하고 선택 후보를 주기적 PBE로 검증할 수 있습니다.
설치와 실제 지원 범위는 [표면 계산 안내](docs/SURFACE_WORKFLOW_KO.md)를 참고하세요.
네 MLIP 전체를 능가했다는 검증 결과는 없습니다. 같은 MACE 모델의 GPU 연산 가속과 독립 DFT 검증을 연결한 도구입니다.

원소·원자 좌표·전체 전하·스핀을 입력해 고립 원자와 분자의 Kohn–Sham DFT를 계산합니다.
기존 수소 전용 계산에서 확장했으며, 물질 이름은 계산에 사용하지 않습니다.
물·H₂·He·OH 등은 편집 가능한 시작 예제입니다.

## 실행

Windows 화면과 WSL DFT 엔진이 설치되어 있으면 **`start_windows.bat`를 더블클릭**하세요.
새 PC에서는 [첫 실행 안내](START_HERE_KO.md)에 따라 Windows Python과 WSL 엔진을 준비하세요.

1. 원소와 x y z 좌표 입력 또는 **XYZ 파일 불러오기**.
2. 좌표 단위, 전체 전하, 스핀 **2S = Nα − Nβ** 확인.
3. **계산 방법** 탭에서 기저함수, ECP, 교환·상관 함수, RKS/UKS 선택.
4. **계산 시작**(F5) → 총에너지·SCF 수렴 여부·밀도 단면 확인. **오비탈** 탭에서 전체 MO 표 확인.
5. 실행 기록에서 설정을 다시 불러오거나 결과 ZIP 내보내기.

전하 0은 중성, +1은 전자 한 개 제거, −1은 전자 한 개 추가입니다.
스핀 입력은 다중도가 아닙니다. singlet은 0, doublet은 1, triplet은 2입니다.
전자 수와 2S의 홀짝이 맞지 않으면 계산을 시작하지 않습니다.
XYZ는 단일 프레임 형식(원자 수, 주석 한 행, 좌표)이며 Å 단위로 읽습니다.
주석의 전하·스핀은 자동 적용하지 않으므로 화면에서 지정하세요.

## 분자 화면의 지원 범위

| 항목 | 현재 지원 |
| --- | --- |
| 계 | 고립 원자·분자·이온, 고정 핵 단일점 계산 |
| 원소 | PySCF H~Og (Z=1~118), 원자 1~100개/800 AO 이하. 자체 H~Kr/160 AO 이하 및 메모리 검사 |
| 전자 상태 | 닫힌껍질 RKS, 열린껍질 UKS; auto는 2S로 결정 |
| SCF 반복 | 자체 Pulay DIIS·밀도 혼합 또는 2차 Newton; 외부 PySCF DIIS/Newton |
| 코어 처리 | 전체전자 또는 무거운 원소의 def2 기저에 대응하는 scalar def2-ECP |
| 기저함수 | STO-3G, 6-31G, 6-31G*, 6-31+G*, cc-pVDZ, def2-SVP/TZVP, Dyall-v2z/v3z |
| 상대론·핵 | Dyall의 auto: 일전자 scalar X2C + Gaussian 유한 핵. 명시적 비상대론/점핵 선택 가능 |
| 교환·상관 | LDA/VWN, PBE, B3LYP, PBE0 |
| 시각화 | 구조 3D 산점도, 총밀도·스핀밀도·HOMO/LUMO의 xy/xz/yz 단면 |
| 결과 | 총에너지, 전자 수, SCF 이력, MO 에너지/점유수/계수, 밀도행렬, CUBE |
| 이 분자 화면의 미지원 | 임의 ECP, 스핀-궤도 결합, IRC, quasi-RRHO·용매 자유에너지, 들뜬 상태; 주기계는 별도 화면 사용 |

UKS와 2S=0을 함께 선택할 수 있지만 깨진 대칭 초기 추정과 여러 스핀 상태 탐색은 자동 제공하지 않습니다.
란타넘족·악티늄족·초중원소를 포함한 118개 원소의 고립 분자·이온·착물을 입력할 수 있습니다.
기저 데이터와 계산 경로의 지원이며, 모든 원소·화합물의 바닥상태나 정량 정확도를 검증했다는 뜻은 아닙니다.
`engines/periodic.py`는 별도 PySCF PBC 엔진이며, 주기계에서는 설치된 GTH 기저/의사퍼텐셜 데이터를 검사합니다.
주기계 GTH MOLOPT-SR 조합은 71개 원소이며 Ag·Au·Pd·Pt를 포함합니다.

| 기저함수 | 설치된 PySCF 2.14.0에서 확인한 지원 원소 |
| --- | --- |
| STO-3G | H~I, 전체전자 최소 기저 |
| 6-31G | H~Zn |
| 6-31G* | H~Kr |
| 6-31+G* | H~Ca |
| cc-pVDZ | K를 제외한 H~Kr |
| def2-SVP, def2-TZVP | H~Rn 중 Ce~Lu 제외, 72개. H~Kr 전체전자, 나머지 대응 ECP |
| Dyall-v2z, Dyall-v3z | H~Og 전체 118개, 전체전자 |

혼합 분자는 포함된 모든 원소를 지원하는 기저만 사용할 수 있습니다. 화면은 가능한 선택지를
갱신하며, 이미 선택한 기저가 맞지 않으면 오류를 표시합니다. 기저를 자동으로 바꾸지 않습니다.

**원소별 기저 (0.8.1):** 화면의 `원소별 기저` 입력란, CLI `--basis-element O=def2-tzvp`(반복 가능),
설정 JSON `"basis_overrides": {"O": "def2-tzvp"}`로 원소마다 다른 기저를 지정합니다. 지정하지 않은 원소는
위 `기저함수`를 사용합니다. 각 원소의 실제 기저 데이터를 검사하며, def2-ECP는 def2 기저가 지정된 무거운 원소에만
붙습니다. 하나라도 Dyall이면 auto 상대론·핵·격자는 X2C/Gaussian/universal이며, ECP 원소와 X2C의 혼합은 차단합니다.
결과에는 `basis_by_element`가 저장되고, 반응·흡착 분석은 **같은 원소의 기저·ECP가 모든 항에서 같을 때만** 계산합니다.
기본값은 def2-SVP이며, 이 기저에 없는 Ce~Lu/Fr~Og는 Dyall 기저를 선택하세요.
이는 이 앱이 검사한 PySCF 내장 기저 데이터의 지원 범위이며, 외부 기저 데이터를 자동으로 받지 않습니다.

## 자체 2차 SCF와 주기계 밴드·DOS (0.9.1)

- **자체 엔진 `scf_solver=newton`**: 궤도 회전 신뢰영역 Newton과 전처리 Steihaug CG. 0.9.1부터 자체 J/K/XC 선형응답과 Libxc 2차 미분으로 해석적 Hessian-벡터 곱을 계산합니다. 점유는 궤도 번호로 유지하며, 실수 내부 불안정 방향을 따라 더 낮은 해를 찾습니다(안정성 허용오차 −1e-4 Ha, 전역 바닥상태 보장 아님). aufbau 위반이면 경고합니다. 기록된 Cr·Co 원자(PBE/STO-3G, grid3, 2S=6/3) 미수렴이 해결됐으며 다른 조건의 수렴은 별도 확인해야 합니다. 물의 같은 궤도에서 응답 계산은 이전 FD 대비 약3.74배 빠릅니다(전체 SCF 배속 아님).
- **주기계 DOS**: 모든 주기계 결과에 k 메쉬 고유값의 가우시안 확장 DOS(`dos_sigma_ev`, 기본 0.1 eV, states/eV/셀)와 HOMO/LUMO·메쉬 갭이 저장됩니다(`dos.csv`).
- **주기계 밴드**: `band_path`(`auto` 또는 `GXWKGL` 같은 고대칭점; CLI `--band-path auto --band-points 60`). 수렴한 밀도로 경로 위 고유값을 비자기일관 계산하고 VBM/CBM·간접/직접 갭을 기록합니다(`bands.csv`). 메쉬 k점에서 밴드값과 SCF 고유값의 차이(Si에서 ~1 meV)를 `mesh_consistency_max_ev`로 남깁니다. graphene slab의 2D GMKG 경로도 실계산 검증했습니다(kz=0, 메쉬 일관성 차3.8e-9 eV).
- **주기계 힘·응력(0.9.2, bulk 전용)**: `integral_scheme=fftdf`(CLI `--integral-scheme fftdf`, 선택 `--fft-mesh NX NY NZ`)와 `calculation=gradient`(CLI `--calculation gradient`). Coulomb과 XC를 같은 균일 FFT 격자에서 계산하며 PySCF 해석적 기울기·응력을 씁니다. 결과 `gradient`: 힘(eV/Å), 응력(Ha/bohr³·GPa, σ = (1/V)∂X/∂ε), 압력(−tr σ/3, GPa), 미분 대상(smearing이면 자유에너지 F = E − TS). FFTDF 에너지는 GDF 에너지와 다르므로 반응 분석에서 섞지 않습니다. slab은 GDF만 가능하므로 힘·응력이 없습니다. FFT mesh·기저·k점 수렴은 직접 확인하세요.

## 자체 SCF 안정성 분석과 수렴 개선 (0.8.6)

- `native_stability`(화면 `수치 설정 > 자체 SCF 안정성 분석`, CLI `--stability check|follow`, 자체 엔진 전용): 수렴한 해가 궤도 회전에 대해 최소인지 검사합니다. 불안정하면 경고하고, `follow`이면 더 낮은 해를 찾아 다시 수렴시킵니다(최대 3회). 실수 내부 회전만 검사합니다(RKS→UKS 삼중항 불안정은 미검사). 결과 `scf_stability`에 회차별 최저 고유값·에너지를 기록합니다.
- 자체 엔진 `auto` 초기 추정은 원자 밀도 중첩(PySCF minao 보조)이며 밀도 혼합은 잔차가 클 때만 적용합니다(`native_damping_residual`, 기본 1e-3). DIIS 반복 수가 약 30~50% 줄었고 Fe 원자가 수렴합니다. Cr·Co의 기록된 정수 점유 사례는 0.9.1 자체 Newton으로 수렴합니다.

## 조화 진동수·ZPE (0.8.5)

`calculation=frequencies`(화면 `계산 종류`, CLI `--calculation frequencies`)는 입력 구조에서 조화 진동수를 계산합니다.
**최적화된 구조에서 실행하세요**(입력 구조의 최대 기울기가 4.5e-4 Ha/bohr를 넘으면 경고).

- Hessian: 해석적 힘의 중심 유한차분(변위 `freq_step_bohr`, 기본 0.005 bohr). 원자 N개이면 기울기 계산 6N+1회이며 각 변위 계산은 `displacements/`에 기록됩니다. 대칭화 전 비대칭 크기를 기록합니다.
- 질량: 최빈 동위원소의 정확한 질량(¹H 1.007825, ¹⁶O 15.994915 …). 병진·회전(선형 분자는 5개) 공간을 정확히 제거해 3N−6(5)개 진동수를 얻습니다. 허수 진동은 음수로 표시합니다.
- 판정: |ν| > 20 cm⁻¹인 허수 진동이 없으면 `minimum`, 1개면 `transition_state_candidate`. 작은 허수 진동은 수치 잡음일 수 있다고 경고합니다.
- ZPE와 진동 열역학(조화 근사: U_vib, S_vib, Cv_vib, 온도 `thermo_temperature_k`, 기본 298.15 K). 비조화성·보정 인자 없음.
- **자유에너지 (0.8.7)**: 이상기체·강체회전자·조화진동자(RRHO)로 병진·회전·전자(다중도 2S+1) 항을 더해 `zpe_corrected_energy_hartree`, `enthalpy_hartree`, `gibbs_free_energy_hartree`를 결과 최상위에 저장합니다. 압력 `thermo_pressure_pa`(기본 101325 Pa), 회전 대칭수 `symmetry_number`(0이면 PySCF 점군 판정으로 자동). 질량은 모든 항에서 최빈 동위원소. 들뜬 전자 상태·quasi-RRHO·용매·표준상태 변환은 없습니다(100 cm⁻¹ 미만 모드가 있으면 경고).
- 반응·흡착 분석(`energy_analysis.py --energy-field gibbs_free_energy_hartree` 등)은 ZPE/H/G를 선택할 수 있고, 이때 모든 항의 온도·압력이 같아야 합니다.
- 격자 응답 항이 없는 힘을 쓰므로 병진·회전 잔여 진동수(물/def2-SVP/격자 3에서 ~17 cm⁻¹)가 생기며, 이 성분은 정확히 제거되고 진단값으로 기록됩니다.

## NEB 최소 에너지 경로 (0.8.9)

`calculation=neb`, `neb_product`(생성물 좌표, 반응물과 **같은 원자 순서**; 화면 `수치 설정`, CLI `--neb-product product.xyz`), `neb_images`(양 끝 포함, 기본 9), `neb_maxg_ev_per_angstrom`(기본 0.05, 평균 기준은 절반).

- geomeTRIC CI-NEB(climbing image)가 사슬을 최적화하고, 이미지마다 에너지·힘은 선택한 엔진이 계산합니다. 양 끝은 고정되므로 **먼저 각각 최적화한 구조**를 넣으세요.
- 초기 경로는 반응물→생성물 직교좌표 선형 보간이며 자동 정렬·원자 재배열을 하지 않습니다. 두 구조의 방향과 원자 순서를 맞춰 주세요.
- 결과: `neb.csv`(이미지별 상대 에너지 kJ/mol), `neb_path.xyz`, `neb/calc_NNNN`(모든 기울기 계산), 상단 폴더에는 최고 에너지 이미지(전이상태 추정)의 전체 결과. `transition_state`로 정밀화하고 진동수로 허수 1개를 확인하는 순서를 권장합니다. IRC는 아직 없습니다.

## 제약·이완 스캔·전이상태 (0.8.8)

- **제약 최적화**: `opt_constraints`(화면 `수치 설정 > 구조 제약`, CLI `--constraints`)에 geomeTRIC 제약 문법을 씁니다. `;`로 줄을 나눌 수 있습니다. 예: `$freeze;xyz 1`(원자 1 고정), `$set;angle 2 1 3 100`(각도 100°로 고정), `$set;distance 1 2 1.5`. 원자 번호는 1부터, 단위는 Å·도입니다. 고정은 geomeTRIC 수렴 허용오차(~1e-4 Å) 안에서 지켜집니다.
- **이완 스캔**: `calculation=scan`, `scan_spec="distance 1 2 0.6 1.0 5"`(또는 `angle i j k …`, `dihedral i j k l …`; 시작·끝·점수). 점마다 해당 좌표를 고정한 최적화를 수행하고 다음 점은 이전 점의 최적 구조에서 시작합니다. `scan/point_NNN/`, `scan.csv`(상대 에너지 kJ/mol 포함), `scan.xyz`. 이완 스캔 최대값은 장벽의 근사이며 전이상태가 아닙니다.
- **전이상태 탐색**: `calculation=transition_state`. geomeTRIC P-RFO가 1차 안장점을 찾으며 시작 Hessian은 우리 엔진 힘의 유한차분(6N회, 결과에 `kind=hessian`으로 기록)입니다. **반드시 진동수 계산으로 허수 진동이 정확히 1개인지 확인하세요.** IRC는 아직 없습니다(NEB는 0.8.9).

## 구조 최적화 (geomeTRIC, 0.8.4)

`calculation=optimize`(화면 `계산 종류`, CLI `--calculation optimize`)는 입력 구조에서 출발해 국소 최소 구조를 찾습니다.
**geomeTRIC**(TRIC 내부좌표, Wang & Song, J. Chem. Phys. 144, 214108, 2016)이 다음 구조를 제안하고,
매 단계의 에너지와 힘은 선택한 엔진(자체/PySCF, D4 포함)이 아래 해석적 기울기로 계산합니다.

- 선택 설치: `wsl --exec .venv-dft/bin/python -m pip install -r requirements-optimize.txt` (검증: geomeTRIC 1.1.1).
- 수렴 기준 `opt_convergence`: `default`(geomeTRIC 기본: ΔE 1e-6 Ha, RMS/최대 기울기 3e-4/4.5e-4 Ha/bohr, RMS/최대 이동 1.2e-3/1.8e-3 Å), `tight`(10배 엄격). 최대 단계 `opt_max_steps`(기본 100).
- 실행 폴더: `steps/step_NNN`(단계별 전체 계산 기록), `optimization.xyz`(궤적), `optimized.xyz`, `final_settings.json`(최종 구조 설정), `geometric.log`. 최종 구조는 요청한 출력 설정으로 한 번 더 계산해 `results.json`에 저장하며 `optimization` 블록에 단계별 에너지·힘과 기준을 기록합니다. `settings.json`은 요청 그대로입니다.
- 단계의 SCF가 수렴하지 않으면 최적화를 중단합니다. 최대 단계 안에 기준을 못 맞추면 상태는 `not_converged`입니다.
- 격자 응답 항이 없는 힘을 쓰므로 `grid_level` 3 이상을 권장합니다(2 이하이면 경고). 결과는 **국소 최소점 후보**이며 진동수 분석은 아직 없습니다. 체크포인트 재시작은 지원하지 않습니다.

## 원자 힘 (해석적 기울기, 0.8.3)

`calculation=gradient`(화면 `계산 종류`, CLI `--calculation gradient`)는 에너지와 함께 입력 구조에서의
원자 힘을 계산합니다. 구조를 바꾸지는 않습니다(구조 최적화는 아직 없음).

- **자체 엔진**: `engines/native_gradient.py`가 해석적 기울기를 직접 조립합니다. 1전자(기저 이동 + 핵 이동), Coulomb·정확교환(2전자 미분 적분을 셸 단위로 계산해 메모리 제한), LDA/GGA XC(AO 2차 도함수), 에너지가중 밀도(Pulay 항), 핵반발로 구성되며, 결과에 성분별로 저장합니다. 미분 적분 값은 PySCF/libcint에서 가져옵니다.
- **PySCF 엔진**: `nuc_grad_method()` (ECP, X2C, Newton 경로 포함).
- D4를 켜면 dftd4의 분산 기울기를 더합니다.
- **격자 응답 항은 포함하지 않습니다**(PySCF 기본값과 같음). 고정 격자 에너지의 정확한 미분이며, 격자가 원자와 함께 움직이는 실제 에너지와의 차이는 물/PBE 기준 격자 수준 1~4에서 6.0e-4, 4.4e-5, 1.1e-5, 1.1e-6 Ha/bohr입니다. 결과의 합력(net force)이 같은 크기의 지표입니다.
- SCF가 수렴하지 않으면 힘을 계산하지 않습니다. 결과 `gradient`: Ha/bohr 기울기, eV/Å 힘(= −기울기), 최대·RMS 힘, 합력, 성분.

## 분산 보정 (DFT-D4, 0.8.2)

`dispersion=d4`(화면 `분산 보정`, CLI `--dispersion d4`)는 Grimme 그룹의 공식 구현 **dftd4**로
DFT-D4 에너지를 계산해 SCF 총에너지에 더합니다. D4는 EEQ 부분전하에 따라 C6가 변하는 모델로,
BJ 감쇠 2체 항과 Axilrod–Teller–Muto 3체 항을 포함합니다(Caldeweyher et al., J. Chem. Phys. 150, 154122, 2019).

- 선택 설치: `wsl --exec .venv-dft/bin/python -m pip install -r requirements-dispersion.txt` (검증: dftd4 4.3.1).
- 공개 감쇠 매개변수가 있는 **PBE / PBE0 / B3LYP**만 허용합니다. LDA에는 D4 매개변수가 없어 차단합니다.
- 지원 원소는 dftd4 참조 데이터가 있는 **H~Lr, Cn~Og (110개)**입니다. Rf~Rg는 차단합니다.
- 고립 분자에서는 라이브러리 기본 실공간 cutoff(2체 60 bohr, 3체 40 bohr)를 **분자 최대 원자간 거리 이상으로 넓혀 모든 쌍을 포함**합니다. 적용 값은 결과에 기록됩니다.
- 분산은 밀도와 무관한 기하 보정입니다. 결과의 `scf_energy_hartree`, `energy_components_hartree.dispersion`, `dispersion`(매개변수·EEQ 전하·버전)에 분리 저장하며, `scf_history`에는 SCF 에너지만 있습니다.
- 반응·흡착 분석은 분산 모델·매개변수·라이브러리 버전이 모든 항에서 같아야 합니다. 주기계 엔진에는 아직 없습니다.

## ECP와 전자 수

`ecp=auto`가 기본값입니다. 무거운 원소에 def2-SVP/TZVP를 쓰면 대응하는 def2-ECP를 적용합니다.
원자별 코어 전자 수는 **Rb~Xe 28개, Cs·Ba·La 46개, Hf~Rn 60개**입니다. 다른 지원 조합은 전체전자 계산입니다.
이는 기저/ECP의 정해진 짝을 적용하는 규칙이며 화면의 전자 수와 ECP 안내에 표시됩니다.
`ecp=none`은 전체전자 기저를 선택할 때만 허용하며, `ecp=def2`는 해당 무거운 원소와 def2 기저가 필요합니다.
무거운 원소의 def2 기저를 ECP 없이 실행하거나 ECP와 전체전자 X2C를 중복 적용하는 설정은 차단합니다.

예를 들어 Zr⁴⁺는 전체 전자 **36개**, ECP 코어 **28개**, 명시적으로 계산하는 전자 **8개**입니다.
RbCl은 전체 54개 중 Rb의 코어 28개를 제외한 26개를 계산합니다.
전하의 정의는 바뀌지 않으며 스핀은 명시적 전자 수와 검사합니다.

ECP 계산의 CUBE·NPZ·밀도 그림은 **코어 전자밀도가 빠진 의사밀도**입니다. 파일 주석과 결과 JSON에도 표시합니다.
원소 번호는 실제 원자번호로 저장하며, 유효 핵전하와 원자별 코어 전자 수는 별도 필드로 기록합니다.
따라서 Zr의 유효 핵전하 12를 원자번호 12인 Mg로 표시하지 않습니다.
ECP와 전체전자 계산의 절대 총에너지를 직접 비교하지 마세요.
Scalar ECP는 코어의 상대론 효과를 반영하지만 이 앱은 스핀-궤도 결합을 계산하지 않습니다.

Dyall 기저는 auto에서 **scalar X2C + Gaussian 유한 핵 + universal 격자**를 사용합니다.
X2C는 일전자 보정이며 SOC와 이전자 picture-change를 포함하지 않습니다. 저장 밀도·오비탈은
X2C 표현으로, 물리적 사성분 밀도에 대한 picture-change 보정이 없습니다.
유한 핵에는 PySCF 기본 질량·경험식을 쓰며 핵 Gaussian 지수를 결과에 기록합니다.
universal 격자는 원소별 반경 표에 의존하지 않는 Gauss-Chebyshev/Becke 적분입니다.
Z>36은 방사 점 수를 PySCF 수준 기본값의 3배로 늘려 실제 점 수를 결과에 기록합니다.
정량 연구에서는 Dyall-v2z/v3z, grid_level, 전자상태를 바꾸어 수렴을 확인해야 합니다.

**범용 입력과 모든 물질의 정확한 계산은 다릅니다.** SCF 수렴은 기저·격자 수렴이나 올바른 전자상태를
보증하지 않습니다. STO-3G는 빠른 검증용이며 음이온에는 확산함수 필요성을 따로 점검하세요.
특히 전이금속은 여러 스핀 상태와 전자배치의 해가 가능하므로 2S를 바꿔 비교해야 합니다.
자동 바닥상태 탐색은 제공하지 않습니다. 지정한 2S, 밀도행렬의 Nα−Nβ, UKS 오비탈의
⟨S²⟩와 목표 S(S+1)의 차이를 결과 JSON에 기록합니다.

외부 PySCF의 DIIS에서 진동하거나 미수렴할 때 **계산 방법 → SCF 반복 방식 → newton**으로 다시 실행할 수 있습니다.
전자 수·2S·XC·기저는 유지하며, 반복 방식은 설정과 결과에 저장합니다. 자동 전환은 하지 않습니다.
Newton의 반복 수는 바깥 반복(macroiteration) 수이며 내부 반복과 구분합니다.
Newton에서 제공되지 않는 밀도 변화 잔차는 JSON의 `null`/CSV 빈칸으로 기록합니다.
Co²⁺ 예제는 Newton을 선택해 둔 예시입니다. 두 반복 방식 모두 미수렴은 성공으로 처리하지 않습니다.

## 결과 해석과 저장

총에너지는 핵 간 반발을 포함합니다. 오비탈 에너지의 합이 아니며,
Kohn–Sham 오비탈 에너지 차이는 실험적 여기 에너지/밴드갭과 동일하지 않습니다.
실행은 `converged`, `not_converged`, `cancelled`, `error`로 구분합니다. 미수렴과 중단은 성공으로 처리하지 않습니다.
분자 계산의 중단 버튼은 안전한 계산 경계에서 반영하며, 저장된 체크포인트가 있으면 실행 기록에서 재시작할 수 있습니다.
미수렴 데이터는 마지막 SCF 반복의 진단값이며 화면·그림·JSON에서 표시합니다.

매 실행은 `outputs/molecular/시각_고유문자/`에 저장됩니다. 기존 실행을 덮어쓰지 않습니다.

| 파일 | 내용 |
| --- | --- |
| settings.json | 좌표·단위·전하·스핀·기저·XC·수치 설정 |
| results.json | 수렴 상태, 총에너지, 오비탈 표, 경고, 단위, 버전 |
| checkpoint.npz | 중간/최종 밀도, 물리 조건 지문과 원본 엔진; 재시작 때 호환성 검사 |
| scf_convergence.png | GUI 실행의 SCF 에너지·잔차 그래프 |
| density.cube | 명시적 계산 전자의 밀도 [electron/bohr³]; ECP에서는 코어 제외 |
| spin_density.cube | UKS의 α−β 스핀밀도 |
| *_homo.cube / *_lumo.cube | 존재하는 각 채널 HOMO/LUMO 진폭 [bohr^(-3/2)] |
| electronic_structure.npz | MO 계수·에너지·점유수, 밀도/AO 겹침 행렬, AO 라벨, 실제 원자번호·유효 핵전하·코어 전자 수 |
| density_grid.npz | 그림용 전자밀도 격자와 상자·원자 좌표 |
| orbitals.csv / scf.csv | 오비탈 표 / SCF 에너지·잔차 이력 |
| molecular.png | 구조와 전자밀도 단면 |
| execution.log | 엔진 실행 로그 |
| error.json | 실패 원인과 설정; 실패 시 생성 |

그림의 좌표는 Å, CUBE 좌표는 bohr입니다. 밀도 단면은 핵 좌표 중심에 가장 가까운 격자면입니다.
SCF의 원자 중심 적분 격자와 CUBE 시각화 격자는 다릅니다. 전자 수는 AO 겹침 행렬을 이용한
`Tr(D S)`로 검사하며, 거친 유한 상자 그림의 적분으로 SCF 정규화를 판단하지 않습니다.
분자 결과 스키마 4(기존 스키마 3과 같은 전자 수 규약)의 `electrons`는 전체 전자 수, `explicit_electrons`는 명시적 계산 전자 수,
`ecp_core_electrons`는 ECP 코어 전자 수입니다. `electron_count_from_density`는 명시적 전자 수와 비교합니다.
`density.representation`은 `all_electron_density`, `ecp_pseudo_density`, `x2c_picture_density`입니다.
`integration_grid`에는 실제 격자 방식·점 수·수치 적분 전자 수 오차를 기록합니다.
`nuclear_repulsion_hartree`는 선택한 모델의 핵 간 반발이며 ECP에서는 유효 핵전하를 사용합니다.
NPZ는 `numpy.load(path, allow_pickle=False)`로 읽을 수 있습니다.

## 명령행

```powershell
.\.venv\Scripts\python.exe main.py --no-show
.\.venv\Scripts\python.exe main.py --config examples\native_water.json --engine native --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\oh.xyz --charge 0 --spin 1 --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\h2.xyz --basis sto-3g --xc pbe --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\water.xyz --basis sto-3g --basis-element O=def2-svp --xc lda,vwn --engine native --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\water.xyz --basis def2-svp --xc pbe --dispersion d4 --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\water_distorted.xyz --basis def2-svp --xc pbe --dispersion d4 --calculation optimize --no-show
.\.venv\Scripts\python.exe main.py --xyz <최적화 결과 폴더>\optimized.xyz --basis def2-svp --xc pbe --dispersion d4 --calculation frequencies --no-show
.\.venv\Scripts\python.exe main.py --config molecule.json --output outputs\study --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\kcl.xyz --charge 0 --spin 0 --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\fe3plus.xyz --charge 3 --spin 5 --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\co2plus.xyz --charge 2 --spin 3 --scf-solver newton --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\rbcl.xyz --charge 0 --spin 0 --ecp auto --no-show
.\.venv\Scripts\python.exe main.py --xyz examples\zr4plus.xyz --charge 4 --spin 0 --ecp auto --no-show
.\.venv\Scripts\python.exe main.py --config examples\auplus.json --no-show
.\.venv\Scripts\python.exe main.py --config examples\og_x2c.json --no-show
```

Windows에서는 `auto`가 프로젝트의 WSL `.venv-dft`를 사용합니다. Linux/macOS는 로컬 PySCF를 사용합니다.
`--backend local|wsl|auto`로 명시할 수도 있습니다. 종료 코드는 수렴 **0**, 실행/입력 실패 **1**, SCF 미수렴 **2**입니다.
`--output`은 실행 폴더들을 저장할 상위 경로입니다. GUI와 CLI는 같은 엔진·입력 검증·출력 형식을 사용합니다.

## 수소 검증기 보존

```powershell
.\.venv\Scripts\python.exe app.py --hydrogen-check
.\.venv\Scripts\python.exe main.py --hydrogen-check --no-show
```

`solver.py`, `config.json`, 기존 방사방정식·해석해 비교와 테스트를 유지했습니다.
수소 검증기는 DFT가 아니며, 중성 He 계산은 새 분자 엔진에서 합니다.
이전 화면과 해석은 [수소 검증기 설명](docs/HYDROGEN_REFERENCE_KO.md),
[방사방정식 참고](MODEL_REFERENCE_KO.md)를 보세요.

## 구현과 검증

`molecular.py`: 입력/검증. `engines/pyscf_engine.py`: 분자 DFT.
`dft_workflow.py`: 로컬/WSL 실행과 기록. `molecular_app.py`: 화면.
`molecular_plotting.py`: 저장된 결과의 시각화.

```powershell
.\.venv\Scripts\python.exe -m unittest -v
wsl --exec .venv-dft/bin/python -m unittest test_molecular -v
$env:QUANTUMSTARTER_GUI_DFT_TEST='1'
.\.venv\Scripts\python.exe -m unittest test_molecular_gui -v
```

PySCF가 없는 Windows Python에서는 실제 물리 테스트가 건너뛰어지므로 WSL 테스트도 실행하세요.
검증: H₂ 에너지, He/H/H₂⁺ 전자 수와 UKS, MO 직교성, 밀도행렬 재구성,
기저/XC 교체, 좌표 단위·평행이동 불변성, 강제 미수렴.
H~Og 확장 검증: 9기저×118원소의 1,062개 조합 중 지원 554개·미지원 508개를 실제 데이터와 대조.
두 Dyall 기저 각각 118개 원소의 분자 구성, 118개 원소의 universal 격자, 모든 지원 def2 ECP 코어를 검사합니다.
선택한 경계 원소의 실제 SCF와 기존 28~40번/전이금속 및 혼합 ECP 분자 RbCl의 계산도 검증합니다.
전자 수·스핀 밀도·코어 전자 수·CUBE/NPZ의 원소 번호·ECP 미수렴도 검사합니다.
이 검증은 모든 금속 착물의 수렴과 정확도를 보증하지 않습니다.
GUI 검증: 입력·XYZ·설정 저장, 실제 WSL 계산, 내보내기, 기록 복원과 오류 복구.

구현 참고: [PySCF DFT](https://pyscf.org/user/dft.html),
[분자/스핀 입력](https://pyscf.org/quickstart.html),
[분자 기저와 ECP](https://pyscf.org/user/gto.html),
[SCF 초기 추정과 안정성](https://pyscf.org/user/scf.html),
[설치 및 Windows/WSL 지원](https://pyscf.org/user/install.html),
[CUBE 출력](https://pyscf.org/_modules/pyscf/tools/cubegen.html).
