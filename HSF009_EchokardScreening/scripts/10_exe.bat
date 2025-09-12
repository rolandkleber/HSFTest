@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"

set "ROOT=%~dp0.."
set "DATA=%ROOT%\data"
set "RESULTS=%ROOT%\results"
set "PY=python"

if not exist "%DATA%" (
  echo ERROR: Data folder not found: %DATA%
  pause & exit /b 1
)
if not exist "%RESULTS%" mkdir "%RESULTS%"

echo === EF from real CSVs (10_ef.py) ===
for /f "delims=" %%F in ('dir /s /b "%DATA%\lv_trace.csv"') do (
  set "CSV=%%~fF"
  for %%S in ("%%~dpF\.") do set "SERIES=%%~nxS"
  echo EF: !SERIES!
  "%PY%" "%ROOT%\scripts\10_ef.py" ^
    --in-csv "!CSV!" ^
    --out "%RESULTS%\!SERIES!" ^
    --fps-list "60,30,20,15,10,8,5"
  echo.
)

echo Done. Check "%RESULTS%\*".
pause
