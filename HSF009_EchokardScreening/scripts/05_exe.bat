@echo off
setlocal EnableExtensions EnableDelayedExpansion
cd /d "%~dp0"

set "ROOT=%~dp0.."
set "RAW=%ROOT%\data\raw"
set "OUT=%ROOT%\data"
set "PY=python"

if not exist "%RAW%" (
  echo ERROR: not found: %RAW%
  pause & exit /b 1
)

echo [EXTRACT] src=%RAW%  dest=%OUT%
"%PY%" "%~dp0\05_extractdata.py" --raw "%RAW%" --out "%OUT%"
echo Done. CSVs in %OUT%\<series>\lv_trace.csv
pause
