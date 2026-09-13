@echo off
rem =====================================================================
rem  Compile the LaTeX tutorial in this folder and copy the PDF to the
rem  Desktop.  This file is intentionally ASCII only so that it can be
rem  parsed correctly regardless of the console code page.
rem  Simply double click this file.
rem =====================================================================
setlocal enabledelayedexpansion
cd /d "%~dp0"

echo.
echo ==========================================================
echo   Compiling the tutorial with xelatex
echo ==========================================================
echo.

where xelatex >nul 2>nul
if errorlevel 1 (
  echo [ERROR] xelatex was not found in PATH.
  echo         Please install MiKTeX or TeX Live first,
  echo         then run this file again.
  echo.
  pause
  exit /b 1
)

set TEXFILE=
for %%f in (*.tex) do set TEXFILE=%%f
if "!TEXFILE!"=="" (
  echo [ERROR] No .tex file found in this folder.
  echo         Current folder: %CD%
  echo.
  pause
  exit /b 1
)

echo Source file : !TEXFILE!
echo Working dir : %CD%
echo.
echo --- pass 1 of 2 ---
xelatex -interaction=nonstopmode "!TEXFILE!"
echo.
echo --- pass 2 of 2 ---
xelatex -interaction=nonstopmode "!TEXFILE!"
echo.

set PDFFILE=
for %%f in (*.pdf) do set PDFFILE=%%f
if "!PDFFILE!"=="" (
  echo [ERROR] No PDF was produced. Compilation failed.
  echo         Open the .log file in this folder and look for lines
  echo         that start with an exclamation mark.
  echo.
  pause
  exit /b 1
)

echo PDF produced : !PDFFILE!
copy /Y "!PDFFILE!" "%USERPROFILE%\Desktop\" >nul
if errorlevel 1 (
  echo [WARN] Could not copy the PDF to the Desktop automatically.
  echo        Please copy it manually.
) else (
  echo The PDF has been copied to your Desktop.
)
echo.
echo All done.
echo.
pause
exit /b 0
