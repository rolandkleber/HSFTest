@echo off
REM === Hard-coded run order, but skip missing/incomplete stages ===

if exist 10_exe.bat (
  call 10_exe.bat || echo [10] skipped (non-critical)
) else (
  echo [10] missing -> skipping
)

if exist 30_exe.bat (
  call 30_exe.bat || (echo [30] FAILED & exit /b 1)
) else (
  echo [30] missing -> skipping & this is IMPORTANT
)

if exist 35_exe.bat (
  call 35_exe.bat || (echo [35] FAILED & exit /b 1)
) else (
  echo [35] missing -> skipping (optional)
)

if exist 40_exe.bat (
  call 40_exe.bat || echo [40] skipped (non-critical)
) else (
  echo [40] missing -> skipping
)

echo [ALL] Done (alpha). See results\
exit /b 0
