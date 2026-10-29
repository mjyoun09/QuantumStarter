@echo off
setlocal
cd /d "%~dp0"
title QuantumStarter - Molecular DFT Engine Setup
wsl --exec python3 -c "import sys; assert sys.version_info >= (3, 10)"
if errorlevel 1 goto missing_wsl
if exist ".venv-dft\pyvenv.cfg" goto install
echo Creating the project-local WSL Python environment...
wsl --exec python3 -m venv .venv-dft
if errorlevel 1 goto failed
:install
echo Installing PySCF and the DFT dependencies in .venv-dft...
wsl --exec .venv-dft/bin/python -m pip install -r requirements-dft.txt
if errorlevel 1 goto failed
wsl --exec .venv-dft/bin/python -c "import pyscf; print('PySCF ready:', pyscf.__version__)"
if errorlevel 1 goto failed
echo Ready. Open start_windows.bat to calculate molecules.
pause
exit /b 0
:missing_wsl
echo A configured WSL Ubuntu installation with Python is required.
echo See START_HERE_KO.md. This script does not change Windows features.
pause
exit /b 1
:failed
echo DFT engine setup failed. Read the error above.
echo If python3-venv is missing, install it inside Ubuntu first.
echo See START_HERE_KO.md.
pause
exit /b 1
