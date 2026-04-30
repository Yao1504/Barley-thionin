@echo off
setlocal enabledelayedexpansion

REM ==================================================
REM Script: run_thionin_pipeline.bat
REM Description:
REM   Build a BLAST protein database, run BLASTP searches
REM   using thionin protein sequences, and filter candidate
REM   barley homologs.
REM
REM Required software:
REM   1. NCBI BLAST+
REM   2. Python >= 3.8
REM   3. Python packages: pandas, biopython
REM
REM Usage:
REM   run_thionin_pipeline.bat ^
REM     INPUT_DIR ^
REM     OUTPUT_DIR ^
REM     SCRIPT_DIR ^
REM     BARLEY_FASTA ^
REM     ARABIDOPSIS_FASTA ^
REM     RICE_FASTA ^
REM     OTHER_FASTA ^
REM     DB_PREFIX ^
REM     OUTPUT_PREFIX
REM
REM Example:
REM   run_thionin_pipeline.bat ^
REM     D:\Thionin\morex\input ^
REM     D:\Thionin\morex\output ^
REM     D:\Thionin\morex\code ^
REM     Hv_Morex.pgsb.Jul2020.aa.fa ^
REM     arabidopsis_4.fasta ^
REM     rice_39.fasta ^
REM     "thionins from 1kb_140.fasta" ^
REM     Hv_Morex_db ^
REM     Hv_thionin
REM ==================================================

REM ==================================================
REM User-defined arguments
REM ==================================================

set INPUT_DIR=%~1
set OUTPUT_DIR=%~2
set SCRIPT_DIR=%~3
set BARLEY_FASTA=%~4
set ARABIDOPSIS_FASTA=%~5
set RICE_FASTA=%~6
set OTHER_FASTA=%~7
set DB_PREFIX=%~8
set OUTPUT_PREFIX=%~9

REM ==================================================
REM Check arguments
REM ==================================================

if "%INPUT_DIR%"=="" (
    echo Error: INPUT_DIR is required.
    goto usage
)

if "%OUTPUT_DIR%"=="" (
    echo Error: OUTPUT_DIR is required.
    goto usage
)

if "%SCRIPT_DIR%"=="" (
    echo Error: SCRIPT_DIR is required.
    goto usage
)

if "%BARLEY_FASTA%"=="" (
    echo Error: BARLEY_FASTA is required.
    goto usage
)

if "%ARABIDOPSIS_FASTA%"=="" (
    echo Error: ARABIDOPSIS_FASTA is required.
    goto usage
)

if "%RICE_FASTA%"=="" (
    echo Error: RICE_FASTA is required.
    goto usage
)

if "%OTHER_FASTA%"=="" (
    echo Error: OTHER_FASTA is required.
    goto usage
)

if "%DB_PREFIX%"=="" (
    echo Error: DB_PREFIX is required.
    goto usage
)

if "%OUTPUT_PREFIX%"=="" (
    echo Error: OUTPUT_PREFIX is required.
    goto usage
)

REM ==================================================
REM Define full paths
REM ==================================================

set BARLEY_FASTA_PATH=%INPUT_DIR%\%BARLEY_FASTA%
set ARABIDOPSIS_FASTA_PATH=%INPUT_DIR%\%ARABIDOPSIS_FASTA%
set RICE_FASTA_PATH=%INPUT_DIR%\%RICE_FASTA%
set OTHER_FASTA_PATH=%INPUT_DIR%\%OTHER_FASTA%
set DB_NAME=%OUTPUT_DIR%\%DB_PREFIX%
set FILTER_SCRIPT=%SCRIPT_DIR%\filter_blast_hits.py

set ARABIDOPSIS_OUT=%OUTPUT_DIR%\arabidopsis_vs_morex.out
set RICE_OUT=%OUTPUT_DIR%\rice_vs_morex.out
set OTHER_OUT=%OUTPUT_DIR%\others_vs_morex.out

REM ==================================================
REM Create output directory
REM ==================================================

if not exist "%OUTPUT_DIR%" (
    mkdir "%OUTPUT_DIR%"
)

REM ==================================================
REM Check input files
REM ==================================================

if not exist "%BARLEY_FASTA_PATH%" (
    echo Error: Barley FASTA not found: %BARLEY_FASTA_PATH%
    exit /b 1
)

if not exist "%ARABIDOPSIS_FASTA_PATH%" (
    echo Error: Arabidopsis FASTA not found: %ARABIDOPSIS_FASTA_PATH%
    exit /b 1
)

if not exist "%RICE_FASTA_PATH%" (
    echo Error: Rice FASTA not found: %RICE_FASTA_PATH%
    exit /b 1
)

if not exist "%OTHER_FASTA_PATH%" (
    echo Error: Other species FASTA not found: %OTHER_FASTA_PATH%
    exit /b 1
)

if not exist "%FILTER_SCRIPT%" (
    echo Error: Python filtering script not found: %FILTER_SCRIPT%
    exit /b 1
)

REM ==================================================
REM Step 1: Build BLAST database
REM ==================================================

echo.
echo ==================================================
echo Step 1: Building BLAST protein database
echo ==================================================

makeblastdb ^
    -in "%BARLEY_FASTA_PATH%" ^
    -dbtype prot ^
    -out "%DB_NAME%"

if errorlevel 1 (
    echo Error: makeblastdb failed.
    exit /b 1
)

REM ==================================================
REM Step 2: Run BLASTP
REM ==================================================

echo.
echo ==================================================
echo Step 2: Running BLASTP searches
echo ==================================================

echo Running BLASTP: Arabidopsis thionins vs barley proteome
blastp ^
    -query "%ARABIDOPSIS_FASTA_PATH%" ^
    -db "%DB_NAME%" ^
    -out "%ARABIDOPSIS_OUT%" ^
    -outfmt 6 ^
    -evalue 1e-5

if errorlevel 1 (
    echo Error: BLASTP failed for Arabidopsis query.
    exit /b 1
)

echo Running BLASTP: Rice thionins vs barley proteome
blastp ^
    -query "%RICE_FASTA_PATH%" ^
    -db "%DB_NAME%" ^
    -out "%RICE_OUT%" ^
    -outfmt 6 ^
    -evalue 1e-5

if errorlevel 1 (
    echo Error: BLASTP failed for Rice query.
    exit /b 1
)

echo Running BLASTP: Other species thionins vs barley proteome
blastp ^
    -query "%OTHER_FASTA_PATH%" ^
    -db "%DB_NAME%" ^
    -out "%OTHER_OUT%" ^
    -outfmt 6 ^
    -evalue 1e-5

if errorlevel 1 (
    echo Error: BLASTP failed for Other species query.
    exit /b 1
)

REM ==================================================
REM Step 3: Filter BLAST hits and extract sequences
REM ==================================================

echo.
echo ==================================================
echo Step 3: Filtering BLAST hits and extracting homologs
echo ==================================================

python "%FILTER_SCRIPT%" ^
    --blast_dir "%OUTPUT_DIR%" ^
    --target_fasta "%BARLEY_FASTA_PATH%" ^
    --output_dir "%OUTPUT_DIR%" ^
    --prefix "%OUTPUT_PREFIX%" ^
    --evalue 1e-5 ^
    --blast_files ^
        "arabidopsis_vs_morex.out" ^
        "rice_vs_morex.out" ^
        "others_vs_morex.out"

if errorlevel 1 (
    echo Error: Python filtering script failed.
    exit /b 1
)

REM ==================================================
REM Finish
REM ==================================================

echo.
echo ==================================================
echo Pipeline finished successfully.
echo Output directory:
echo %OUTPUT_DIR%
echo ==================================================

exit /b 0

:usage
echo.
echo Usage:
echo   run_thionin_pipeline.bat INPUT_DIR OUTPUT_DIR SCRIPT_DIR BARLEY_FASTA ARABIDOPSIS_FASTA RICE_FASTA OTHER_FASTA DB_PREFIX OUTPUT_PREFIX
echo.
echo Example:
echo   run_thionin_pipeline.bat D:\Thionin\morex\input D:\Thionin\morex\output D:\Thionin\morex\code Hv_Morex.pgsb.Jul2020.aa.fa arabidopsis_4.fasta rice_39.fasta "thionins from 1kb_140.fasta" Hv_Morex_db Hv_thionin
echo.
exit /b 1
