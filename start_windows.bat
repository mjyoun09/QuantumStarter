@echo off
setlocal
cd /d "%~dp0"
title QuantumStarter - Setup and Launch
if exist ".venv\Scripts\python.exe" goto environment_ready

echo Creating a local Python environment...
if exist "%LocalAppData%\Programs\Python\Python313\python.exe" (
    "%LocalAppData%\Programs\Python\Python313\python.exe" -m venv .venv
    goto environment_created
)
py -3.13 -c "import sys" >nul 2>&1
if not errorlevel 1 (
    py -3.13 -m venv .venv
    goto environment_created
)
py -3.12 -c "import sys" >nul 2>&1
if not errorlevel 1 (
    py -3.12 -m venv .venv
    goto environment_created
)
py -3.11 -c "import sys" >nul 2>&1
if not errorlevel 1 (
    py -3.11 -m venv .venv
    goto environment_created
)
python -c "import sys; assert sys.version_info >= (3, 11)" >nul 2>&1
if errorlevel 1 goto missing_python
python -m venv .venv

:environment_created
if not exist ".venv\Scripts\python.exe" goto failed

:environment_ready
".venv\Scripts\python.exe" -c "import numpy, scipy, matplotlib, ase" >nul 2>&1
if not errorlevel 1 goto check_tk
echo Installing the desktop dependencies, including ASE structure readers.
echo This step requires an internet connection.
".venv\Scripts\python.exe" -m pip install -r requirements.txt
if errorlevel 1 goto failed

:check_tk
".venv\Scripts\python.exe" -c "import tkinter" >nul 2>&1
if errorlevel 1 goto missing_tk
echo Opening QuantumStarter...
".venv\Scripts\python.exe" app.py %*
if errorlevel 1 goto failed
exit /b 0

:missing_python
echo Python 3.11 or newer is required.
echo Install Python with Tcl/Tk support, then run this file again.
echo See START_HERE_KO.md for Korean instructions.
pause
exit /b 1

:missing_tk
echo Tkinter is missing. Modify your Python installation to include Tcl/Tk.
echo See START_HERE_KO.md for Korean instructions.
pause
exit /b 1

:failed
echo Setup or launch failed. The error message is shown above.
echo See START_HERE_KO.md for Korean troubleshooting instructions.
pause
exit /b 1
