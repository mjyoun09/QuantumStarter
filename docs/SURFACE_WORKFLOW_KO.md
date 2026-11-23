# 고체·표면: 빠른 탐색 + DFT 검증

0.6은 고립 분자와 별도로 PySCF PBC 엔진을 실행합니다. 이전에 고체·표면이 안 된 이유는 격자, 주기 경계,
k점, 금속 점유수와 적절한 의사퍼텐셜 입력이 없었기 때문입니다. 입력 제한 문구만 없애서는 고칠 수 없는 부분입니다.

첨부된 네 노트북은 MACE-MP, SevenNet, CHGNet, MatterSim이라는 학습된 원자간 포텐셜을 사용합니다.
전자밀도를 매번 푸는 DFT와 목적·비용이 다릅니다. 현재 구현은 새 포텐셜을 학습한 것이 아니라,
MACE로 후보를 고르고 독립적인 주기적 DFT로 선택 구조를 검증하는 워크플로입니다.
네 모델 전체보다 우수하다는 성능 주장은 하지 않습니다.

## 화면에서 실행

1. `start_windows.bat` → **고체·표면 / GPU 탐색**, 또는 `start_windows.bat --periodic`.
2. **주기적 DFT** 탭에서 Si 벌크/Pt(111) 소형 예제를 선택하거나 POSCAR/CIF/extxyz를 불러옵니다.
3. `bulk`는 xyz 주기, `slab`은 xy 주기 및 z 비주기 2D Coulomb 경계입니다. slab의 k점은 `kx ky 1`입니다.
   POSCAR/CIF에는 표면 경계 정보가 없으므로, 불러오기 전 `slab`을 선택해 명시하세요.
4. 격자 벡터는 행별 입력이며 좌표와 함께 Å 단위입니다. 스핀은 **입력 셀당** `Nα−Nβ`입니다.
5. **계산 시작**. E, F, (E+F)/2, SCF 수렴 여부와 원시 결과 경로가 표시됩니다.

0.9.3부터 bulk / FFTDF의 원자 힘·응력 및 원자/셀 완화를 GUI에서 선택할 수 있습니다. **Si 셀 완화**, FFT mesh 검사, 중단 요청은 [bulk 셀 계산 안내](PERIODIC_CELL_KO.md)를 참조하세요. 아래 표면 후보 검증은 계속 `gdf` / `single_point` 설정을 사용합니다.

Si/Pt 예제의 낮은 격자 수준과 작은 kmesh는 동작 확인용입니다. 연구용 정밀도가 아닙니다.
일반적인 시작값은 코드 기본값(PBE / DZVP-MOLOPT-SR-GTH, 적분 격자 3)이지만 이것도 정확도 보증은 아닙니다.

## 표면 H 후보 탐색

별도의 `.venv-mlip` 환경을 사용합니다. Windows에서는 `setup_mlip_windows.bat`로 설치합니다.
CUDA 가능한 NVIDIA GPU가 필요하며 CLI는 `--device cpu`도 허용합니다. 모델을 받는 최초 실행에는 인터넷이 필요합니다.
MACE 모델은 프로젝트 `.model_cache/mace`에 저장합니다. GPU가 없거나 모델 로딩에 실패해도 다른 모델로 바꾸지 않습니다.

**표면 H 탐색 → DFT 검증** 탭에서:

1. HEA 또는 단일 금속, x/y 크기와 층 수, seed를 지정합니다.
2. **후보 생성 + GPU 탐색**: 각 표면 격자 위치의 ontop/bridge/fcc/hcp에서 H를 시작해 완화합니다.
   표면 원자는 고정하며 H만 움직입니다. 기본 3×3×4는 표면 36원자 + H이고 후보는 36개입니다.
   HEA seed 0의 조성은 Ag8 Au7 Cu7 Pd7 Pt7입니다. 첨부 노트북의 실제 생성 코드와 같습니다.
3. H의 최종 위치가 주기적 최소 이미지 거리 0.15 Å 이내인 후보를 묶고, 그룹마다 최저 에너지 구조를 남깁니다.
   서로 다른 후보 중 에너지가 낮은 3개와 나머지 중 seed로 고른 탐색 후보 1개를 선택합니다.
4. **선택 후보 DFT 검증**: 주기적 DFT 탭의 기저·적분·SCF 설정과 xy k점을 사용하고 z k점은 1입니다.
   GUI의 표면 검증은 원자가 전자의 홀짝에 맞춘 최소 스핀(짝수 0/홀수 1)을 사용하고 결과에 기록합니다.
   바닥상태 스핀 탐색이 아니며, 다른 스핀은 CLI 설정 파일로 명시해 비교하세요.
5. 선택된 동일 조성·셀의 구조끼리 상대 에너지와 순위를 비교합니다. 미수렴/오류 결과는 비교에서 제외합니다.

`initial_site`는 시작 위치입니다. H가 다른 사이트로 이동하거나 여러 출발점이 같은 위치에 도달할 수 있습니다.
고정 표면, 특정 합금 배치 및 coverage에 대한 탐색입니다. 전체 표면 배치/스핀의 최적해를 찾았다는 뜻은 아닙니다.
현재 검증은 **MLIP 완화 구조에서 DFT 단일점**이며 DFT 구조 최적화가 아닙니다.

## 명령행

Windows PowerShell, 프로젝트 폴더:

```powershell
.\.venv\Scripts\python.exe main.py --config examples/pt111_small.json --no-show
.\.venv\Scripts\python.exe main.py --structure POSCAR --boundary slab --kmesh 3 3 1 --spin 1 --no-show
.\.venv\Scripts\python.exe surface_workflow.py generate --output outputs/surfaces/my_hea --material HEA --size 3 3 4 --seed 0
wsl --exec .venv-mlip/bin/python -u surface_workflow.py screen outputs/surfaces/my_hea/candidates.json --top-k 3 --explore 1 --acceleration cueq
```

탐색이 출력한 `screen_.../screen.json` 경로로 검증합니다:

```powershell
.\.venv\Scripts\python.exe surface_workflow.py verify outputs/surfaces/my_hea/screen_실행시각/screen.json --backend wsl
```

`--settings settings.json`으로 `kmesh`, `spin`, `basis`, `grid_level`, `precision`, `conv_tol`, `smearing_ev`, `max_cycle`,
`threads`, `memory_mb` 등을 지정할 수 있습니다. `examples/periodic_smoke_settings.json`은 저비용 기능 검사 전용입니다.
`--limit 2`로 선택 후보 중 2개만 계산할 수 있습니다. 전체 37원자 HEA의 주기적 DFT는 MLIP보다 훨씬 무거우며
GPU는 MACE에 사용하고 현재 PySCF PBC 엔진은 CPU에서 실행합니다.

`--model sevennet|chgnet|mattersim`은 선택 어댑터입니다. 별도 호환 환경에 해당 패키지가 필요합니다.
기본 설치와 실제 실행 검증 대상은 **MACE-MP medium / float64**입니다. 나머지 어댑터의 실행·성능은 이 릴리스에서 검증하지 않았습니다.
SevenNet은 7net-0, CHGNet은 0.3.0, MatterSim은 v1.0.0-5M으로 고정하며 실패 시 다른 모델로 대체하지 않습니다.

## 측정한 GPU 가속

화면의 cuEquivariance 옵션은 같은 MACE-MP medium 가중치와 float64를 유지하고 GPU 연산 구현을 바꿉니다.
RTX 5080 Laptop GPU에서 37원자 HEA+H 구조 3개를 각각 7회 측정한 에너지+힘 평가의 중앙값은
기본 e3nn보다 **2.60배** 빨랐습니다. 최대 에너지 차이는 2.85×10⁻¹⁴ eV 미만,
힘 성분 차이는 1.34×10⁻¹⁵ eV/Å 미만이었습니다.
이는 워밍업 이후 같은 모델의 추론 비교입니다. 설치/모델 로딩/완화/DFT를 포함한 전체 처리속도나
다른 세 모델보다 정확하다는 뜻은 아닙니다. 원시 측정은 `outputs/mace_acceleration_benchmark.json`입니다.
공식 구현 설명: [MACE CUDA 가속](https://mace-docs.readthedocs.io/en/latest/guide/cuda_acceleration.html).

재측정:

```powershell
wsl --exec .venv-mlip/bin/python benchmark_mace.py outputs/surfaces/my_hea/ontop_0_0.extxyz outputs/surfaces/my_hea/fcc_1_1.extxyz outputs/surfaces/my_hea/hcp_2_2.extxyz --output outputs/mace_acceleration_benchmark_new.json --repeats 7
```

## 물리적 지원 범위와 결과

- 주기계: 중성 셀, 3D 벌크 또는 xy 평면 슬래브, PBE/GTH-PBE, KRKS/KUKS, Fermi smearing.
- 기저: GTH DZVP-MOLOPT-SR, SZV-MOLOPT-SR, DZVP, TZVP. 모든 기저에 모든 원소가 있는 것은 아닙니다.
  실제 설치 데이터가 없는 조합은 명확히 거부합니다. 목표 금속의 MOLOPT-SR/PBE 조합을 검사했습니다.
- 실제 원자가 전자: Cu/Ag/Au 11, Pd/Pt 18, H 1. 이 값을 전체 원자번호와 구분합니다.
- 원자 200개, AO 2000개, k점 216개는 입력 상한이지 이 PC에서 실용적인 규모를 보장하는 값이 아닙니다.
- bulk FFTDF는 원자 힘·응력과 셀/원자 동시 최적화를 지원합니다. 이 표면 워크플로의 slab/GDF 경로에는 힘·응력·구조 최적화가 없습니다. 주기계 SOC·분산·하이브리드 XC·하전 표면은 후속 과제입니다. [bulk 지원 조건](PERIODIC_CELL_KO.md)
- SCF가 마지막 추가 검사에서 실패하면 남은 반복 한도 안에서 같은 조건·최종 밀도로 한 번 재시도합니다.
  허용 오차를 느슨하게 만들지 않으며 그래도 실패하면 `not_converged`입니다.
- 각 실행의 `settings.json`, `results.json`, `scf.csv`, `scf.chk`, `electronic_structure.npz`, `execution.log`를 보관합니다.
  에너지는 입력 셀당 값이며 k점 수로 곱하지 않습니다. `E`, `F=E−TS`, `(E+F)/2`를 별도로 저장합니다.
- `screen.json`은 모델 버전/장치/가중치 지문, 각 구조 지문, 원자 힘의 벡터 노름, 완화 수렴과 시간을 기록합니다.
  DFT 검증 전에 구조 변경·조성/셀 불일치를 검사합니다.
- 비교는 DFT 최저 후보를 공통 기준으로 잡은 상대 에너지 MAE와 0.01 eV보다 큰 차이의 순서 일치율입니다.
  DFT 선택 표본 밖의 정확도를 뜻하지 않으며 원래 모델끼리 절대 총에너지 오프셋을 비교하지 않습니다.

절대 수소 흡착에너지는 `E(slab+H) − E(slab) − ½E(H₂)`의 각 항을 일관된 계산 모델로 얻어야 합니다.
현재 사이트 비교에는 H₂ 기준을 넣지 않았습니다. 고립 H₂의 MLIP 오차가 사이트 비교에 섞이는 것을 피하면서,
절대 흡착에너지인 것처럼 표시하지 않습니다. 정량 연구에는 k점·기저·진공·슬래브 두께·smearing 수렴,
적절한 스핀/자성 상태와 물리적 보정의 검토가 필요합니다.

## 구현 근거

[PySCF 결정 구조·2D 경계·k점·GTH](https://pyscf.org/user/pbcgto.html),
[PySCF 주기 SCF·smearing·자유에너지](https://pyscf.org/user/pbc/scf.html),
[MACE 공식 저장소](https://github.com/ACEsuit/mace).
