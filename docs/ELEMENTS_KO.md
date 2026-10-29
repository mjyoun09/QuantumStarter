# 118개 원소의 입력과 실제 계산 범위 · 0.7

분자 화면은 **H부터 Og까지 118개 원소**를 받습니다. 상단 **원소별 지원표**에서 각 원소에 맞는 기저와 주기계 지원을 확인할 수 있습니다. 좌표·전하·스핀을 직접 지정하며, 지원 데이터가 없는 조합은 계산 전에 거부합니다.

## 분자 기저와 ECP

| 기저 | 원소 수 | 범위 |
| --- | ---: | --- |
| STO-3G | 53 | H~I |
| 6-31G | 30 | H~Zn |
| 6-31G* | 36 | H~Kr |
| 6-31+G* | 20 | H~Ca |
| cc-pVDZ | 35 | H~Kr, K 제외 |
| def2-SVP / def2-TZVP | 각각 72 | H~Rn, Ce~Lu 제외 |
| Dyall-v2z / Dyall-v3z | 각각 118 | H~Og 전체전자 |

표는 설치된 PySCF 2.14.0의 실제 데이터에서 생성했습니다. 이름만 인식하도록 확장한 것이 아닙니다.
def2에서 Rb~Xe 코어28, Cs·Ba·La 코어46, Hf~Rn 코어60개를 대체합니다. Au⁺는 전체78전자 중 18개를 명시적으로 계산합니다. Dyall은 코어를 포함한 전체전자 기저입니다.

`data/element_capabilities.json`은 Windows 화면에서도 읽을 수 있습니다. 재생성 명령:

```powershell
wsl --exec .venv-dft/bin/python scripts/build_element_catalog.py
```

갱신 후 실제 데이터 대조 테스트를 수행하세요. 설치된 라이브러리의 기저가 바뀌어도 지원 범위를 추측하지 않습니다.

## 상대론, 핵 모형, 적분

- `relativity=auto`: Dyall에서 `sfx2c1e()` 일전자 scalar X2C. 다른 기저에서는 기존 동작 유지. 명시적 `none`/`x2c` 선택 가능.
- `nuclear_model=auto`: Dyall에서 Gaussian 유한 핵. PySCF 기본 동위원소 질량 및 경험식으로 생성하며 실제 Gaussian 지수 `nuclear_gaussian_zeta_per_bohr2`를 저장합니다. 사용자 동위원소 지정은 아직 없습니다.
- `grid_scheme=auto`: Dyall에서 universal. Gauss-Chebyshev 방사 격자, 각도 축소 없는 Lebedev 격자, 반경 보정 없는 Becke 분할입니다. Z>36의 방사 점 수는 해당 PySCF 수준 기본값의 3배로 늘립니다(예: Og level2 270점, level3 315점). 실제 원소별 점 수는 `integration_grid.atom_grid_overrides`에 기록합니다. 원소 반경 값을 임의로 추가하지 않습니다. 기존 기저는 Treutler 방식을 유지합니다.
- PySCF 기본 Treutler 방사 매개변수는 Lr까지만 있습니다. Rf~Og에서 이를 강제로 선택하면 사전 오류를 표시합니다.

X2C에 **SOC와 이전자 picture-change는 포함되지 않습니다.** 밀도·오비탈은 X2C 표현이며 물리적 사성분 밀도로의 picture-change는 수행하지 않습니다. 결과 JSON, CUBE 주석, 그림에 표현을 구분합니다. ECP와 전체전자 X2C를 동시에 켜는 설정은 차단합니다.

`Tr(DS)` 전자 수 외에 SCF 수치 격자에서 적분한 전자 수와 오차도 기록합니다. 격자 전자 수 오차가 0.01개보다 크면 경고합니다. 이 기준을 통과해도 에너지의 격자 수렴을 보증하지 않습니다. `grid_level`을 높여 에너지를 직접 비교해야 합니다.

근거: [PySCF Dyall 기저](https://pyscf.org/pyscf_api_docs/pyscf.gto.basis.dyall-basis.html), [scalar 상대론](https://pyscf.org/user/scf.html#scalar-relativistic-correction), [방사 격자 구현](https://pyscf.org/_modules/pyscf/dft/radi.html), [격자 및 분할 구현](https://pyscf.org/_modules/pyscf/dft/gen_grid.html). 이 앱이 선택한 auto 규칙과 검사 기준은 앱의 설계이며 문헌 정확도 보증이 아닙니다.

## 주기계 및 ML

| 주기계 PBE/GTH-PBE 조합 | 지원 원소 수 |
| --- | ---: |
| GTH DZVP-MOLOPT-SR / SZV-MOLOPT-SR | 각각 71 |
| GTH DZVP | 22 |
| GTH TZVP | 18 |

분자 Dyall의 118개 원소 지원은 주기계 의사퍼텐셜이나 ML 체크포인트의 지원 원소를 늘리지 않습니다. 주기계는 실제 기저와 GTH-PBE 양쪽 데이터가 있는 원소만 허용합니다. ML 모델은 체크포인트별 학습 범위와 검증이 필요합니다. 없는 원소를 다른 원소로 바꾸거나 미학습 결과를 검증된 예측으로 표시하지 않습니다.

## 실행과 검증

```powershell
.\.venv\Scripts\python.exe main.py --config examples\auplus.json --no-show
.\.venv\Scripts\python.exe main.py --config examples\og_x2c.json --no-show
wsl --exec .venv-dft/bin/python -m unittest test_molecular test_periodic -v
wsl --exec .venv-dft/bin/python scripts/validate_heavy_elements.py
```

지원표는 9기저×118원소, 1,062조합을 대조하며 실제 지원 554조합을 확인합니다. Dyall 두 기저 각각 118원소의 분자 구성과 universal 격자 생성도 검사합니다. **118종 모두의 SCF·바닥상태·정량 정확도를 검증한 것은 아닙니다.** 실제 경계 원소 계산과 조건·결과는 [검증 기록](VALIDATION_KO.md)에 별도로 기록합니다.

분자 최대 100원자/800 AO, 지정 전하·스핀, 고정 구조 단일점이라는 기존 자원·방법 제한은 유지합니다. 초중원소·악티늄족의 연구 결과에는 SOC, 전자상관, 기저·격자, 전자배치에 대한 추가 검토가 필요합니다.
