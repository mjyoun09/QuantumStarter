@echo off
setlocal
cd /d "%~dp0"
if exist ".venv-mlip\pyvenv.cfg" goto install
wsl --exec python3 -m venv .venv-mlip
if errorlevel 1 goto failed
:install
wsl --exec .venv-mlip/bin/python -m pip install torch==2.7.1 --index-url https://download.pytorch.org/whl/cu128
if errorlevel 1 goto failed
wsl --exec .venv-mlip/bin/python -m pip install -r requirements-mlip.txt
if errorlevel 1 goto failed
wsl --exec .venv-mlip/bin/python -c "import torch; print('CUDA:', torch.cuda.is_available()); import mace"
if errorlevel 1 goto failed
echo Ready. Open start_windows.bat and select the surface workflow.
pause
exit /b 0
:failed
echo MLIP setup failed. Read the error above. See docs/SURFACE_WORKFLOW_KO.md.
pause
exit /b 1
