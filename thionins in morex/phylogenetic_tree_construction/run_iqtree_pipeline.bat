@echo off
setlocal enabledelayedexpansion

REM ============================================================
REM Script: run_iqtree_pipeline.bat
REM
REM Description:
REM   Automated IQ-TREE phylogenetic analysis pipeline for
REM   protein sequence alignments.
REM
REM   The pipeline performs:
REM     1. ModelFinder model selection
REM     2. Maximum-likelihood tree construction
REM     3. 1000 bootstrap replicates
REM     4. Generation of a tree with branch support cutoff >= 50%
REM     5. Logging of analysis progress
REM
REM Required software:
REM   - IQ-TREE 3
REM
REM Required input:
REM   - Trimmed protein alignment in FASTA format
REM
REM Output:
REM   - <PREFIX>_MF.iqtree
REM   - <PREFIX>.treefile
REM   - <PREFIX>_cut50.treefile
REM   - <PREFIX>.iqtree
REM   - <PREFIX>.ckp.gz
REM   - <PREFIX>_bootstrap_progress.log
REM
REM Usage:
REM   run_iqtree_pipeline.bat INPUT_ALIGNMENT OUTPUT_DIR PREFIX THREADS
REM
REM Example:
REM   run_iqtree_pipeline.bat ^
REM     input/modified_protein_align_trimmed.fasta ^
REM     output_tree ^
REM     thionin_ML ^
REM     6
REM ============================================================


REM ============================================================
REM User-defined arguments
REM ============================================================

set "INPUT=%~1"
set "OUTPUT_DIR=%~2"
set "PREFIX_NAME=%~3"
set "THREADS=%~4"


REM ============================================================
REM Default settings
REM ============================================================

if "%THREADS%"=="" (
    set "THREADS=4"
)

if "%PREFIX_NAME%"=="" (
    set "PREFIX_NAME=thionin_ML"
)


REM ============================================================
REM Check arguments
REM ============================================================

if "%INPUT%"=="" goto usage
if "%OUTPUT_DIR%"=="" goto usage

if not exist "%INPUT%" (
    echo Error: input alignment file not found:
    echo %INPUT%
    exit /b 1
)


REM ============================================================
REM Check IQ-TREE executable
REM ============================================================

where iqtree3 >nul 2>nul

if errorlevel 1 (
    echo Error: iqtree3 was not found in PATH.
    echo Please install IQ-TREE 3 and add it to your system PATH.
    exit /b 1
)


REM ============================================================
REM Create output directory
REM ============================================================

if not exist "%OUTPUT_DIR%" (
    mkdir "%OUTPUT_DIR%"
)


REM ============================================================
REM Define output prefix and log file
REM ============================================================

set "PREFIX=%OUTPUT_DIR%\%PREFIX_NAME%"
set "LOGFILE=%OUTPUT_DIR%\%PREFIX_NAME%_bootstrap_progress.log"


REM ============================================================
REM Initialize log file
REM ============================================================

echo ========================================================== > "%LOGFILE%"
echo IQ-TREE Phylogenetic Analysis Pipeline >> "%LOGFILE%"
echo Date: %date%  Time: %time% >> "%LOGFILE%"
echo Input alignment: %INPUT% >> "%LOGFILE%"
echo Output prefix: %PREFIX% >> "%LOGFILE%"
echo Threads: %THREADS% >> "%LOGFILE%"
echo ========================================================== >> "%LOGFILE%"


REM ============================================================
REM Step 1: ModelFinder
REM ============================================================

if not exist "%PREFIX%_MF.iqtree" (
    echo [%time%] Running ModelFinder... >> "%LOGFILE%"

    iqtree3 ^
        -s "%INPUT%" ^
        -st AA ^
        -m MFP ^
        -nt %THREADS% ^
        -pre "%PREFIX%_MF" >> "%LOGFILE%" 2>&1

    if errorlevel 1 (
        echo Error: ModelFinder failed. See log file:
        echo %LOGFILE%
        exit /b 1
    )
) else (
    echo [%time%] ModelFinder result already exists: %PREFIX%_MF.iqtree >> "%LOGFILE%"
)


REM ============================================================
REM Step 2: Extract best-fit model
REM ============================================================

set "BEST_MODEL="

for /f "tokens=6" %%a in ('findstr /C:"Best-fit model according to BIC:" "%PREFIX%_MF.iqtree"') do (
    set "BEST_MODEL=%%a"
)

if "%BEST_MODEL%"=="" (
    echo Error: could not detect best-fit model from:
    echo %PREFIX%_MF.iqtree
    echo Please check the ModelFinder output file.
    exit /b 1
)

echo [%time%] Best model detected: %BEST_MODEL% >> "%LOGFILE%"
echo Best model detected: %BEST_MODEL%


REM ============================================================
REM Step 3: ML tree + 1000 bootstrap
REM ============================================================

echo [%time%] Starting ML tree + 1000 bootstrap... >> "%LOGFILE%"

iqtree3 ^
    -s "%INPUT%" ^
    -st AA ^
    -m %BEST_MODEL% ^
    -b 1000 ^
    -T %THREADS% ^
    -pre "%PREFIX%" >> "%LOGFILE%" 2>&1

if errorlevel 1 (
    echo Error: ML tree / bootstrap analysis failed. See log file:
    echo %LOGFILE%
    exit /b 1
)


REM ============================================================
REM Step 4: Generate cutoff >=50%% bootstrap tree
REM ============================================================

if exist "%PREFIX%.treefile" (
    echo [%time%] Generating cutoff >=50%% bootstrap tree... >> "%LOGFILE%"

    iqtree3 ^
        -t "%PREFIX%.treefile" ^
        -minsup 0.5 ^
        -pre "%PREFIX%_cut50" ^
        -redo >> "%LOGFILE%" 2>&1

    if errorlevel 1 (
        echo Error: bootstrap cutoff tree generation failed. See log file:
        echo %LOGFILE%
        exit /b 1
    )
) else (
    echo Error: ML tree file not found:
    echo %PREFIX%.treefile
    exit /b 1
)


REM ============================================================
REM Finish
REM ============================================================

echo [%time%] Analysis complete! >> "%LOGFILE%"
echo Output files: >> "%LOGFILE%"
echo   %PREFIX%.treefile          - ML tree with 1000 bootstrap >> "%LOGFILE%"
echo   %PREFIX%_cut50.treefile    - ML tree with bootstrap >=50%% >> "%LOGFILE%"
echo   %PREFIX%.iqtree            - main analysis report >> "%LOGFILE%"
echo   %PREFIX%_MF.iqtree         - ModelFinder report >> "%LOGFILE%"
echo   %PREFIX%.ckp.gz            - checkpoint file for resume >> "%LOGFILE%"
echo ========================================================== >> "%LOGFILE%"

echo.
echo ==========================================================
echo IQ-TREE analysis finished successfully.
echo Output directory:
echo %OUTPUT_DIR%
echo Log file:
echo %LOGFILE%
echo ==========================================================

pause
exit /b 0


:usage
echo.
echo Usage:
echo   run_iqtree_pipeline.bat INPUT_ALIGNMENT OUTPUT_DIR PREFIX THREADS
echo.
echo Example:
echo   run_iqtree_pipeline.bat input/modified_protein_align_trimmed.fasta output_tree thionin_ML 6
echo.
exit /b 1
