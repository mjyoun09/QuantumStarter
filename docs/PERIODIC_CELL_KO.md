# Bulk 셀 완화와 FFT 격자 수렴 (0.9.3)

`bulk` / `integral_scheme=fftdf` / PBE-GTH에서 원자 위치와 셀을 함께 완화한다. PySCF가 SCF·원자 힘·응력을 계산하고 ASE의 BFGS와 FrechetCellFilter가 최적화한다. 외부 압력은 0 GPa로 고정이다. slab/GDF, 원자 고정 제약, 외부 압력 지정, 최적화 재시작은 지원하지 않는다.

## 화면과 CLI

고체·표면 화면의 **Si 셀 완화** 예제는 `examples/silicon_cell.json`을 불러온다. `계산 종류=cell_optimize`, `적분 방식=fftdf`를 선택한다. `gradient`는 입력 구조의 힘·응력만 계산한다. `cell_relaxation=full`은 셀 모양과 체적, `isotropic`은 모양을 유지한 체적을 완화하며 둘 다 원자를 움직인다. 설정 영역은 세로 스크롤로 이동한다. 밴드 경로는 비워 두면 생략된다.

Windows 프로젝트 폴더에서:

```powershell
.venv\Scripts\python.exe main.py --config examples/silicon_cell.json --no-show
.venv\Scripts\python.exe main.py --config examples/silicon_cell.json --cell-relaxation full --cell-max-steps 40 --no-show
```

WSL에서는 `.venv-dft/bin/python`으로 같은 Python 파일을 실행한다. `--no-show`는 그림 파일을 저장하고 창을 띄우지 않는다.

CLI/JSON 수렴 설정은 `--cell-force-tol` / `cell_force_tol` (기본 0.02 eV/Å), `--cell-stress-tol-gpa` / `cell_stress_tol_gpa` (0.1 GPa), `--cell-energy-tol-ev` / `cell_energy_tol_ev` (1e-5 eV/cell), `--cell-max-steps` / `cell_max_steps` (50)다. 예제는 최대 40단계다.

## 수치 정의와 완료 판정

- smearing이 0이면 총에너지 E, 양수이면 전자 Mermin 자유에너지 F=E−TS를 최소화한다. `(E+F)/2`를 최소화하지 않는다. 전자 smearing의 F는 원자 진동·격자 열역학을 포함한 물질의 전체 자유에너지가 아니다.
- FFT 정수 격자는 초기 구조에서 한 번 정하고 모든 단계에서 고정한다. `fft_mesh`가 비었으면 PySCF가 `precision`으로 정한 초기 mesh를 사용한다. 셀이 커지면 실제 공간 격자 간격도 커진다. k점은 동일한 **상대 좌표** 메쉬이며 Cartesian k점은 역격자에 맞춰 바뀐다.
- 최대 원자 힘은 원자별 힘 벡터의 노름 중 최댓값이다. 응력은 σ=(1/V)∂E/∂ε 또는 smearing 시 ∂F/∂ε이며, 압력은 −tr(σ)/3이다. `full`은 대칭화한 응력 성분의 최대 절댓값, `isotropic`은 |tr(σ)/3|을 판정한다. 등방 완화에서는 전단·편차 응력이 남을 수 있다.
- 힘·응력과 연속 단계 E/F 변화의 세 기준을 모두 검사한다. 초기 구조에는 이전 단계가 없으므로 에너지 변화 조건을 생략하지만 최종 재계산 일치는 검사한다. 최종 구조를 독립 SCF로 다시 계산하여 SCF·힘·응력·에너지를 재검사한다.
- 최대 단계에 도달하면 `not_converged`다. 중간 SCF 미수렴은 힘을 최적화기에 넘기지 않고 오류로 중단한다. GUI **계산 중단 요청**은 SCF 콜백/단계 경계에서 처리하며, 실행 중인 힘·응력 적분이 끝날 때까지 잠시 걸릴 수 있다. `cancelled`는 성공으로 처리하지 않는다.

힘은 PySCF의 해석적 원자 기울기다. PySCF 2.14 `pbc.grad.krks_stress` / `kuks_stress`의 응력은 운동에너지·겹침 적분 일부에 strain 1e-5 중심 유한차분을 사용한다. 전체 응력을 순수 해석적 구현으로 부르지 않는다. 셀 좌표는 [ASE FrechetCellFilter 문서](https://ase.gitlab.io/ase/ase/filters.html)를 따른다.

## FFT mesh 수렴 검사

GUI에서 `21 21 21;25 25 25;29 29 29`처럼 mesh를 3개 이상 입력하고 **FFT mesh 수렴 확인**을 누른다. 각 성분은 감소하면 안 되고 매 단계 적어도 한 성분은 증가해야 한다. 셀 완화 설정을 불러온 경우에도 검사는 **입력 구조를 고정한 gradient 계산**이다.

```powershell
.venv\Scripts\python.exe periodic_convergence.py --config examples/silicon_cell.json --meshes "21 21 21;25 25 25;29 29 29"
```

마지막 **연속 두 구간**의 세 차이가 모두 기준 이하여야 통과한다.

| 비교량 | 기본 허용오차 | CLI 조정 |
| --- | --- | --- |
| 원자당 E/F 변화의 절댓값 | 0.001 eV/atom | `--energy-tol` |
| 원자별 힘 차이 노름의 최댓값 | 0.01 eV/Å | `--force-tol` |
| 응력 성분 차이의 최대 절댓값 | 0.1 GPa | `--stress-tol` |

모든 요청 mesh를 실행하며 마지막 mesh를 `selected_mesh`로 기록한다. GUI의 현재 mesh를 자동 변경하지 않는다. 최종 구조의 `final_settings.json`으로도 검사를 반복한다. 기저와 k점 수렴은 이 검사에 포함되지 않는다. 적은 수의 mesh에서 차이가 작은 것은 수치 점검 결과이며 모든 더 조밀한 mesh의 수렴을 증명하지 않는다.

## 보존되는 결과

- `settings.json`: 최초 요청. `steps/step_NNNN/`: 각 구조의 설정·SCF·힘·응력·밀도.
- `cell_history.json`, `cell_optimization.csv`, `cell_optimization.extxyz`: 단계별 E/F·힘·응력·체적과 구조 이력. JSON/구조 궤적은 매 단계 저장한다.
- `final/`: 최종 구조 독립 재계산. 루트에는 전체 최적화 판정이 끝난 후 `results.json`을 쓴다. 중단 도중 단일점 수렴을 최적화 성공으로 남기지 않는다.
- `final_settings.json`, `optimized.extxyz`: 마지막 구조. **미수렴 실행에서도 생길 수 있으므로 반드시 `results.json.status`를 함께 확인한다.** `final_settings.json`의 계산 종류는 `gradient`다. 다시 완화하려면 `--calculation cell_optimize`를 지정한다.
- `results.json.cell_optimization`: 알고리즘·버전, E/F 정의, 고정 mesh, 기준별 판정. `final_geometry_config`: 실제 최종 에너지의 구조. 그래프와 반응 에너지 호환성 검사도 이 최종 셀을 사용한다.
- mesh 검사는 별도 고유 실행 폴더에 `mesh_convergence.json`/CSV와 모든 `runs/`를 보존한다. 오류·취소·SCF 미수렴은 수렴 선택을 만들지 않는다.

실측과 재현 명령은 [검증 기록](VALIDATION_KO.md)의 0.9.3 절을 참조한다. SZV/2³ k점 Si 예제는 작동 검증용이며 실리콘의 정밀 평형 격자나 실험 일치 자료가 아니다. 임의 구조의 수렴, 전역 최소점, 포논 안정성은 보장하지 않는다.
