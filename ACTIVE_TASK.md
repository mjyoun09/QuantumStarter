# 현재 작업

## 작업 식별 정보

- 작업 ID: QS-20261029-NATIVE-DFT
- 작업명: 첨부 범용 DFT 요구사항의 누락 구현, 자체 SCF 기반부터 단계별 검증
- 담당 AI: ChatGPT Codex (단독 담당)
- 상태: HANDOFF_READY

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

정확한 다음 구현 지점:

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
