# 처음 실행하기 · 분자·주기계 DFT

## 현재 PC에서 실행

**`start_windows.bat`를 더블클릭**하면 범용 분자 DFT 화면이 열립니다.
Windows 화면용 `.venv`와 WSL 계산용 `.venv-dft`를 따로 사용합니다.
기본 예제는 물 분자입니다. 좌표와 전하·스핀을 확인한 뒤 **계산 시작**을 누르세요.
자체 SCF를 확인하려면 시작 예제에서 **자체 DFT · 물 (LDA 검증)**을 선택하세요.
`native`는 자체 밀도·SCF 엔진(H~Kr), `pyscf`는 외부 엔진(기저별 최대118원소)입니다.
**SCF 실시간** 탭에 에너지·밀도 변화·잔차가 나타납니다. 중단 버튼과 실행 기록의 체크포인트 재시작을 사용할 수 있습니다.
구조 탭의 미리보기에서 거리/각도를 확인하고 수치 설정의 프리셋에서 실제 변경 값을 확인하세요.
새 기능의 라이브러리 의존성과 검증 범위는 [자체 엔진 안내](docs/NATIVE_DFT_KO.md)를 참고하세요.
XYZ 파일을 불러온 경우 전하와 스핀은 직접 확인하세요. 주석의 값을 자동 적용하지 않습니다.

**원소 확장:** 분자 계산은 H~Og(118개)를 지원합니다. 상단 **원소별 지원표**를 확인하세요.
기본 def2 기저에 없는 원소는 **계산 방법 → dyall-v2z**를 선택합니다. auto에서 scalar X2C와 유한 핵을 적용합니다.
시작 예제 **오가네손 · Og (118번 · X2C)** 또는 `examples/og_x2c.json`으로 확인할 수 있습니다.
주기계·ML 모델은 원소 범위가 별도입니다. [원소별 설명](docs/ELEMENTS_KO.md)을 참고하세요.

**고체·표면:** 상단 **고체·표면 / GPU 탐색** 버튼 또는 `start_windows.bat --periodic`으로 엽니다.
주기계는 POSCAR/CIF/extxyz, 격자, k점, 셀당 스핀을 사용합니다. Si 벌크와 Pt(111) 작동 확인 예제가 있습니다.
빠른 표면 탐색에는 **`setup_mlip_windows.bat`**로 별도 GPU 모델 환경을 설치합니다.
이 PC에서는 DFT 및 MACE/CUDA 환경을 설치했습니다. 자세한 순서와 제한은 [표면 계산 안내](docs/SURFACE_WORKFLOW_KO.md)를 참고하세요.

## 새 Windows PC에 설치

1. Python 3.11 이상을 Tcl/Tk 포함으로 설치합니다. 3.12 또는 3.13을 사용할 수 있습니다.
2. `start_windows.bat`를 실행합니다. 화면용 패키지는 프로젝트 `.venv`에 설치됩니다.
3. DFT 계산에는 WSL의 Linux 환경이 필요합니다. PySCF는 Windows Python에서 직접 지원되지 않습니다.
   [공식 PySCF 설치 문서](https://pyscf.org/user/install.html)
4. WSL이 없는 PC에서는 관리자 PowerShell에서 아래 명령을 실행하고 필요한 재부팅/Ubuntu 초기 설정을 마칩니다.

```powershell
wsl --install -d Ubuntu-24.04
```

Ubuntu 터미널에서 Python 가상환경 도구를 준비합니다.

```bash
sudo apt update
sudo apt install python3-venv
```

5. 프로젝트 폴더의 **`setup_dft_windows.bat`**를 실행합니다.
   프로젝트 `.venv-dft`에 PySCF를 설치하며 Windows 기능을 자동 변경하지 않습니다.
6. **`start_windows.bat`**로 화면을 열고 계산합니다.
7. (선택) D4 분산 보정을 쓰려면 PowerShell에서 한 번 실행합니다. 없으면 `분산 보정 none` 계산은 그대로 동작합니다.
   `wsl --exec .venv-dft/bin/python -m pip install -r requirements-dispersion.txt`

최초 설치에는 인터넷이 필요하며, 계산 자체는 로컬에서 수행합니다.
한글이나 공백이 포함된 프로젝트 경로도 Windows에서 WSL 경로로 변환해 실행합니다.

## Linux / macOS

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-dft.txt
.venv/bin/python app.py
```

Linux 화면에는 배포판의 `python3-tk` 패키지와 그래픽 환경이 필요할 수 있습니다.
창 없이 실행하려면 `.venv/bin/python main.py --no-show`를 사용하세요.

## 입력 예

물 분자, 좌표 단위 angstrom, 전하 0, 스핀 0:

```text
O 0 0 0
H 0 -0.757 0.587
H 0 0.757 0.587
```

OH 라디칼은 `examples/oh.xyz`를 읽고 전하 0, 스핀 1로 설정합니다.
스핀은 **2S = Nα − Nβ**이며 다중도가 아닙니다.
미수렴은 성공한 계산이 아닙니다. 계산 방법과 수치 설정을 검토해 다시 실행하세요.
현재 H~Og (원자번호 1~118)의 고립계와 해당 원소를 포함한 착물을 지원합니다.
고체·표면은 별도 주기계 화면을 사용하세요. Ag·Au·Pd·Pt는 분자·주기계 양쪽에서 지원합니다.
전이금속은 여러 스핀 상태를 비교해야 하며, 수렴한 상태가 바닥상태라고 자동 판정하지 않습니다.

시작 예제에서 KCl, Fe³⁺(전하 3, 2S=5), Co²⁺(전하 2, 2S=3)를 선택할 수 있습니다.
Co²⁺ 시작 예제는 수렴을 위해 SCF 반복 방식 `newton`을 사용합니다.
Cu⁺, Kr, RbCl, Zr⁴⁺ 외에 Au⁺와 Og 예제를 추가했습니다.
XYZ 파일로 읽으면 전하와 2S는 직접 입력하세요. 기본 def2-SVP는 Ce~Lu를 제외한 H~Rn을 지원합니다.
지원되지 않는 원소·기저 조합은 계산 전에 차단하며, 자세한 표는 README를 참고하세요.

무거운 원소에 def2 기저를 쓰면 **ECP auto**가 원소별 코어 전자 28/46/60개를 대신 처리합니다.
Zr⁴⁺ 예제는 전체 전자 36개 중 8개를 직접 계산합니다. 화면에 두 전자 수를 구분해 표시합니다.
ECP 계산의 밀도 그림과 파일에는 코어 전자밀도가 없습니다. 전체전자 계산과 절대 총에너지를 직접 비교하지 마세요.

## 문제 해결

- **DFT용 WSL 환경 없음**: `setup_dft_windows.bat`를 먼저 실행하세요.
- **WSL 접근 실패**: Ubuntu 초기 설정과 `wsl --list --verbose`를 확인하세요.
- **ensurepip/python3-venv 오류**: Ubuntu 안에서 `sudo apt install python3-venv`를 실행하세요.
- **Windows Python 별칭 오류**: Microsoft Store 별칭 대신 설치된 Python의 실제 경로를 사용하세요.
  시작 스크립트는 `%LocalAppData%\Programs\Python\Python313\python.exe`도 찾습니다.
- **Tkinter 없음**: Python 설치에서 Tcl/Tk를 추가하세요.
- **계산 실패**: 실행 로그 탭 또는 `outputs/molecular/실행폴더/error.json`을 확인하세요.
- **전자 수와 스핀 불일치**: 전하와 2S의 홀짝 조건을 확인하세요.
- **ECP 설정 오류**: 무거운 원소의 def2 기저는 ECP auto 또는 def2가 필요합니다. none으로 끌 수 없습니다.
- **Ce·U·Og 등의 기저 오류**: Dyall-v2z/v3z를 선택하세요. PySCF 2.14 이상이 필요합니다.
- **SCF가 계속 미수렴**: 계산 방법 탭에서 반복 방식을 `newton`으로 바꿔 재시도할 수 있습니다.
  수렴해도 다른 스핀 상태보다 에너지가 낮은지는 따로 비교해야 합니다.
- **음이온 결과가 불안정함**: 확산함수가 포함된 기저와 기저/격자 수렴을 점검하세요.

기존 수소 검증 화면은 `start_windows.bat --hydrogen-check` 또는
`.venv\Scripts\python.exe app.py --hydrogen-check`로 열 수 있습니다.
결과 파일과 과학적 한계는 [README_KO.md](README_KO.md)를 보세요.
