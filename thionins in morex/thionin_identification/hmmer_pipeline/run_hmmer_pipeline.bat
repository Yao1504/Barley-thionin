@echo off
setlocal enabledelayedexpansion

REM ==================================================
REM Script: run_hmmer_pipeline.bat
REM Description:
REM   Run HMMER hmmscan using a thionin HMM profile against
REM   the barley protein dataset, then filter HMMER hits and
REM   generate a README file.
REM
REM Note:
REM   This script is designed for HMMER 3.0 Windows/Cygwin-style builds.
REM   Windows paths are converted to Cygwin-style paths for hmmscan.
REM
REM Requirements:
REM   - HMMER 3.0 for Windows/Cygwin
REM   - Python
REM   - PowerShell
REM
REM Usage:
REM   run_hmmer_pipeline.bat HMMER_DIR INPUT_DIR OUTPUT_DIR SCRIPT_DIR
REM
REM Example:
REM   run_hmmer_pipeline.bat ^
REM     C:\tools\hmmer-3.0-windows ^
REM     input ^
REM     output_hmm ^
REM     scripts
REM ==================================================


REM ==================================================
REM User-defined paths
REM ==================================================

set HMMER_DIR=%~1
set INPUT_DIR=%~2
set OUTPUT_DIR=%~3
set SCRIPT_DIR=%~4


REM ==================================================
REM Default file names
REM Modify these names if your input files are different
REM ==================================================

set HMM_NAME=PF00321.hmm
set PROTEIN_NAME=Hv_Morex.pgsb.Jul2020.aa.fa
set TBL_NAME=PF00321_vs_morex.tbl


REM ==================================================
REM Check arguments
REM ==================================================

if "%HMMER_DIR%"=="" goto usage
if "%INPUT_DIR%"=="" goto usage
if "%OUTPUT_DIR%"=="" goto usage
if "%SCRIPT_DIR%"=="" goto usage


REM ==================================================
REM Define Windows paths
REM ==================================================

set HMM_FILE_WIN=%INPUT_DIR%\%HMM_NAME%
set PROTEIN_FILE_WIN=%INPUT_DIR%\%PROTEIN_NAME%
set TBL_OUT_WIN=%OUTPUT_DIR%\%TBL_NAME%

set FILTER_SCRIPT=%SCRIPT_DIR%\filter_hmm_hits_cygwin.py
set README_SCRIPT=%SCRIPT_DIR%\generate_hmm_readme.ps1


REM ==================================================
REM Create output directory
REM ==================================================

if not exist "%OUTPUT_DIR%" (
    mkdir "%OUTPUT_DIR%"
)


REM ==================================================
REM Check required files
REM ==================================================

if not exist "%HMMER_DIR%\hmmpress.exe" (
    echo Error: hmmpress.exe not found:
    echo %HMMER_DIR%\hmmpress.exe
    exit /b 1
)

if not exist "%HMMER_DIR%\hmmscan.exe" (
    echo Error: hmmscan.exe not found:
    echo %HMMER_DIR%\hmmscan.exe
    exit /b 1
)

if not exist "%HMM_FILE_WIN%" (
    echo Error: HMM profile file not found:
    echo %HMM_FILE_WIN%
    exit /b 1
)

if not exist "%PROTEIN_FILE_WIN%" (
    echo Error: protein FASTA file not found:
    echo %PROTEIN_FILE_WIN%
    exit /b 1
)

if not exist "%FILTER_SCRIPT%" (
    echo Error: Python filtering script not found:
    echo %FILTER_SCRIPT%
    exit /b 1
)

if not exist "%README_SCRIPT%" (
    echo Error: PowerShell README script not found:
    echo %README_SCRIPT%
    exit /b 1
)


REM ==================================================
REM Convert Windows paths to Cygwin-style paths
REM
REM Note:
REM   This conversion assumes paths are on a drive-letter path,
REM   such as C:\project\input\file.hmm.
REM ==================================================

set HMM_FILE=%HMM_FILE_WIN:\=/%
set PROTEIN_FILE=%PROTEIN_FILE_WIN:\=/%
set TBL_OUT=%TBL_OUT_WIN:\=/%

set HMM_FILE=/cygdrive/%HMM_FILE:~0,1%%HMM_FILE:~2%
set PROTEIN_FILE=/cygdrive/%PROTEIN_FILE:~0,1%%PROTEIN_FILE:~2%
set TBL_OUT=/cygdrive/%TBL_OUT:~0,1%%TBL_OUT:~2%


REM ==================================================
REM Step 1: Press HMM file
REM ==================================================

echo.
echo ============================
echo Step 1: Pressing HMM file
echo ============================

"%HMMER_DIR%\hmmpress.exe" "%HMM_FILE%"

if errorlevel 1 (
    echo Error: hmmpress failed.
    exit /b 1
)


REM ==================================================
REM Step 2: Run HMMER hmmscan
REM ==================================================

echo.
echo ============================
echo Step 2: Running HMMER hmmscan
echo ============================

"%HMMER_DIR%\hmmscan.exe" ^
    --tblout "%TBL_OUT%" ^
    "%HMM_FILE%" ^
    "%PROTEIN_FILE%"

if errorlevel 1 (
    echo Error: hmmscan failed.
    exit /b 1
)


REM ==================================================
REM Step 3: Filter HMMER hits
REM ==================================================

echo.
echo ============================
echo Step 3: Filtering hits and generating FASTA / CSV
echo ============================

python "%FILTER_SCRIPT%" ^
    "%TBL_OUT_WIN%" ^
    "%PROTEIN_FILE_WIN%" ^
    "%OUTPUT_DIR%"

if errorlevel 1 (
    echo Error: HMMER filtering script failed.
    exit /b 1
)


REM ==================================================
REM Step 4: Create README
REM ==================================================

echo.
echo ============================
echo Step 4: Creating README file
echo ============================

powershell -ExecutionPolicy Bypass ^
    -File "%README_SCRIPT%" ^
    "%OUTPUT_DIR%"

if errorlevel 1 (
    echo Error: README generation failed.
    exit /b 1
)


REM ==================================================
REM Finish
REM ==================================================

echo.
echo ============================
echo HMM scan finished successfully.
echo Output directory:
echo %OUTPUT_DIR%
echo ============================

pause
exit /b 0


:usage
echo.
echo Usage:
echo   run_hmmer_pipeline.bat HMMER_DIR INPUT_DIR OUTPUT_DIR SCRIPT_DIR
echo.
echo Example:
echo   run_hmmer_pipeline.bat C:\tools\hmmer-3.0-windows input output_hmm scripts
echo.
exit /b 1
