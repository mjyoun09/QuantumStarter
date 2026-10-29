이 문서는 초기 수치 모델과 명령행 사용법을 보존한 참고 문서입니다.
현재 사용 화면과 실행법은 README_KO.md를 먼저 확인하세요.

# QuantumStarter v0.1

수소 원자의 에너지와 방사형 파동함수를 직접 수치 계산하는 Python 시작 프로젝트입니다.
향후 DFT를 구현하기 전, 해밀토니안 구성과 고유값 계산을 확인하는 첫 단계입니다.

**현재 버전은 DFT가 아닙니다.** 전자 간 상호작용, 교환·상관 함수, SCF 반복,
분자 구조, 촉매 표면, 흡착에너지 계산은 아직 들어 있지 않습니다.
에너지를 해석식으로 찍어내는 방식이 아니라, 유한차분 행렬을 만들어 수치적으로 풉니다.
행렬의 고유값·고유벡터 계산에는 SciPy를 사용합니다.

## 1. VS Code에서 실행하기

Python 3.11 또는 3.12와 VS Code의 Microsoft Python 확장을 준비하세요.
VS Code에서 압축 해제한 **QuantumStarter 폴더 자체**를 엽니다.
`main.py`, `solver.py`, `config.json`이 바로 보이는 폴더입니다.

VS Code에서 터미널을 열고 아래 명령을 한 줄씩 실행하세요.
Windows PowerShell 기준이며, 가상환경을 활성화하지 않아도 실행됩니다.

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe main.py
```

Python 3.11만 설치했다면 첫 줄의 `-3.12`를 `-3.11`로 바꾸세요.
`py` 명령이 없지만 `python --version`으로 3.11/3.12가 확인되면
첫 줄 대신 `python -m venv .venv`를 사용하세요.

이후 VS Code에서 `Ctrl+Shift+P` → `Python: Select Interpreter`를 실행하고
`.venv`의 Python을 선택하면 `main.py`의 실행 버튼이나 F5로 실행할 수 있습니다.
F5 설정은 포함되어 있습니다. 디버깅 확장을 요청하면 Microsoft Python Debugger를 설치하세요.
`solver.py`는 계산 모듈이고, 실행할 파일은 `main.py`입니다.

macOS/Linux 터미널에서는 다음과 같습니다.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python main.py
```

첫 실행에서 그래프 글꼴 캐시를 만들며 잠시 기다릴 수 있습니다.
그래프 창을 닫으면 실행이 끝납니다. 창을 띄울 수 없는 환경에서는 다음처럼 실행하세요.

```powershell
.\.venv\Scripts\python.exe main.py --no-show
```

## 2. 현재 출력

기본 설정은 수소 한 원자의 1s, 2s, 3s 상태입니다. 즉 `l=0`으로 고정하고
`n=1,2,3`을 계산합니다. 2p, 3p 등 다른 각운동량 상태는 별도로 설정해야 합니다.
콘솔에 수치 에너지, 해석적 에너지, 절대오차가 표시됩니다.
그래프의 실선은 수치해이고 점선은 같은 상태의 해석해입니다.

`outputs` 폴더에 아래 파일을 저장합니다.

| 파일 | 내용 |
| --- | --- |
| `results.json` | 설정, 모델 가정, 라이브러리 버전, 에너지, 오차, 정규화, 평균 반지름 |
| `energies.csv` | 상태별 에너지와 오차. 엑셀로 열 수 있음 |
| `radial.csv` | 거리별 축약 방사함수 u, 방사확률밀도 P, 해석해 u |
| `wavefunctions.png` | 파동함수와 방사확률밀도 그래프 |

같은 출력 폴더로 재실행하면 네 결과 파일을 갱신합니다. 비교 결과를 남기려면
`--output outputs/grid12000`처럼 다른 경로를 지정하세요.
기본 출력 경로는 터미널의 위치와 관계없이 프로젝트 안의 `outputs`입니다.

## 3. 값을 바꾸는 곳

`config.json`을 수정하고 다시 실행하세요. JSON에는 주석을 넣지 마세요.

| 설정 | 기본값 | 의미 |
| --- | --- | --- |
| `nuclear_charge` | 1 | 핵전하 Z. 전자 수는 항상 한 개 |
| `angular_momentum_l` | 0 | 방위 양자수. 0은 s, 1은 p |
| `num_states` | 3 | 해당 l에서 낮은 순서대로 구할 방사형 상태 수 |
| `r_max_bohr` | 60.0 | 방사형 계산영역의 바깥 경계, 단위 bohr |
| `grid_points` | 6000 | 양 끝 경계를 제외한 내부 격자점 수 |

Z를 2로 바꾸면 전자 한 개인 He+를 계산합니다. 중성 헬륨이 아닙니다.
같은 l에서 첫 번째 상태의 주양자수는 `n=l+1`입니다. 예를 들어 l=1이면 2p부터 시작합니다.
각도 방향의 구면조화함수와 3차원 오비탈 시각화는 아직 구현하지 않았습니다.
6000~12000점 정도로 시작하세요. 격자점을 너무 늘리는 것보다
영역 크기와 격자 간격을 따로 바꾸며 정확도를 확인하는 것이 우선입니다.

## 4. 계산 원리와 단위

원자 단위계에서 다음 방정식을 풉니다.

```text
[-1/2 d²/dr² + l(l+1)/(2r²) - Z/r] u(r) = E u(r)
u(r) = r R(r)
u(0) = u(r_max) = 0
```

핵은 움직이지 않는 무한 질량 점핵으로, 전자는 비상대론적으로 취급합니다.
고정된 l에서 바깥 경계를 유한한 위치에 둔 방사형 수치 문제입니다.
실제 수소의 유한 핵질량 효과나 미세구조는 포함하지 않습니다.

두 번째 미분을 2차 중앙차분으로 바꾸면 삼중대각 해밀토니안이 됩니다.
`solver.py`가 이 행렬을 구성하고 SciPy의 `eigh_tridiagonal`로 낮은 고유상태를 구합니다.
해석식 `E_n = -Z²/(2n²)`는 계산 후 검증에만 사용합니다.

에너지 내부 단위는 hartree이고 거리 단위는 bohr입니다.
1 hartree는 약 27.2114 eV, 1 bohr는 약 0.529177 Å입니다.
정확한 변환 상수는 설치된 SciPy의 물리 상수에서 가져옵니다.

그래프의 `P(r)=|u(r)|²`는 방사확률밀도입니다.
`P(r)dr`는 반지름 r~r+dr인 구껍질에서 전자를 발견할 확률이고,
3차원 부피당 확률밀도 `|psi|²`와 다릅니다.
수소 1s에서 방사확률밀도의 최대는 약 1 bohr, 평균 반지름은 약 1.5 bohr입니다.
수치해는 `sum(|u|²)*dr=1`로 정규화합니다. 경계의 u는 0입니다.
파동함수의 전체 부호는 임의이므로 그림 비교를 위해 해석해와 부호만 맞춥니다.

## 5. 검증 실행

이 프로그램은 수치 계산 코드이므로 이미 알려진 해와 비교하는 검증을 포함합니다.

```powershell
.\.venv\Scripts\python.exe -m unittest -v
```

검증 항목은 수소 1s·2s·3s 에너지, 정규화와 직교성,
1s 파동함수와 평균 반지름, 격자 세분화에 따른 2차 수렴,
He+ 2p 에너지, 유한 영역 오차, 잘못된 입력 처리입니다.
기본 격자에서 1s의 에너지 절대오차는 약 0.00034 eV입니다.
이것은 고정된 무한 질량 핵이라는 동일한 모델의 해석해에 대한 오차이며,
실험값에 대한 정확도나 이후 DFT의 정확도를 뜻하지 않습니다.

격자 수렴 실험은 r_max를 고정하고 grid_points만 6000 → 12000으로 늘려 비교하세요.
경계 수렴 실험은 r_max와 grid_points를 함께 두 배로 해서 간격을 거의 유지하고 비교하세요.
들뜬 상태는 더 큰 계산영역이 필요합니다. 에너지가 양수인 상태를 결합상태로 해석하지 마세요.

## 6. 다음에 이어갈 위치

`solver.py`는 물리 계산, `plotting.py`는 그래프,
`main.py`는 설정 읽기·실행·결과 저장을 담당합니다.
검증된 수소 계산기를 유지하면서 다음 단계는 별도 `helium_scf.py`에 구현하는 것을 권합니다.
먼저 구대칭 헬륨의 전자밀도와 Hartree 퍼텐셜을 다루고,
교환·상관 근사와 SCF 수렴 조건을 추가해야 실제 Kohn–Sham DFT 단계가 됩니다.
현재 수소 계산만으로 촉매 성능을 예측했다고 해석할 수는 없습니다.

VS Code의 AI에 이어서 설명할 때 다음 내용을 사용해도 됩니다.

> 이 프로젝트는 DFT 엔진을 직접 공부하며 만들기 위한 시작 코드다.
> 현재는 한 전자 원자의 방사형 슈뢰딩거 방정식을 2차 유한차분으로 푼다.
> 우선 README_KO.md, solver.py, test_solver.py를 읽고 현 구조와 단위를 설명해라.
> 기존 수소 계산과 검증을 유지하며, 다음 단계인 구대칭 헬륨 SCF/DFT의
> 물리적 가정과 구현 순서를 먼저 정리해라. 촉매나 분자 기능을 바로 덧붙이지 마라.

## 참고 자료

- [SciPy 공식 문서: eigh_tridiagonal](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh_tridiagonal.html)
- [University of Texas: 수소 원자의 에너지와 방사형 방정식](https://farside.ph.utexas.edu/teaching/qm/lectures/node42.html)
- [OpenStax: The Hydrogen Atom](https://openstax.org/books/university-physics-volume-3/pages/8-1-the-hydrogen-atom)

작성 환경에서 Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0, Matplotlib 3.10.8로 검증했습니다.
Windows의 그래프 창과 VS Code F5 동작은 사용자의 설치 환경에 따라 달라질 수 있으며,
작성 환경에서는 Linux에서 창 없이 계산·파일 출력 및 그래프 저장을 검증합니다.
