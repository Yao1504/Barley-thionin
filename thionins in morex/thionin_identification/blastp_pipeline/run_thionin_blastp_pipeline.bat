@echo off
setlocal enabledelayedexpansion

REM ==================================================
REM Script: run_thionin_pipeline.bat
REM Description:
REM   Run BLASTP searches and downstream filtering to
REM   identify thionin homologs in barley reference line Morex.
REM
REM Requirements:
REM   - NCBI BLAST+
REM   - Python (with pandas, biopython)
REM   - PowerShell
REM
REM Usage:
REM   run_thionin_pipeline.bat INPUT_DIR OUTPUT_DIR SCRIPT_DIR
REM
REM Example:
REM   run_thionin_pipeline.bat input output scripts
REM ==================================================


REM ==================================================
REM User-defined paths (FROM COMMAND LINE)
REM ==================================================
set INPUT_DIR=%~1
set OUTPUT_DIR=%~2
set SCRIPT_DIR=%~3

REM Default file names (can be modified by user if needed)
set BARLEY_FASTA=Hv_Morex.pgsb.Jul2020.aa.fa
set ARABIDOPSIS_FASTA=arabidopsis_4.fasta
set RICE_FASTA=rice_39.fasta
set OTHER_FASTA=thionins_from_1kb_140.fasta

set DB_NAME=%OUTPUT_DIR%\Hv_Morex_db


REM ==================================================
REM Check arguments
REM ==================================================
if "%INPUT_DIR%"=="" goto usage
if "%OUTPUT_DIR%"=="" goto usage
if "%SCRIPT_DIR%"=="" goto usage


REM ==================================================
REM Create output folder if not exists
REM ==================================================
if not exist "%OUTPUT_DIR%" (
    mkdir "%OUTPUT_DIR%"
)


REM ==================================================
REM Step 1: Build BLAST database
REM ==================================================
echo.
echo ==================================================
echo Step 1: Building BLAST database
echo ==================================================

makeblastdb ^
    -in "%INPUT_DIR%\%BARLEY_FASTA%" ^
    -dbtype prot ^
    -out "%DB_NAME%"


REM ==================================================
REM Step 2: Run BLAST
REM ==================================================
echo.
echo ==================================================
echo Step 2: Running BLASTP searches
echo ==================================================

echo Running BLAST for Arabidopsis thionins...
blastp ^
    -query "%INPUT_DIR%\%ARABIDOPSIS_FASTA%" ^
    -db "%DB_NAME%" ^
    -out "%OUTPUT_DIR%\arabidopsis_vs_morex.out" ^
    -outfmt 6 ^
    -evalue 1e-5


echo Running BLAST for Rice thionins...
blastp ^
    -query "%INPUT_DIR%\%RICE_FASTA%" ^
    -db "%DB_NAME%" ^
    -out "%OUTPUT_DIR%\rice_vs_morex.out" ^
    -outfmt 6 ^
    -evalue 1e-5


echo Running BLAST for Other thionins...
blastp ^
    -query "%INPUT_DIR%\%OTHER_FASTA%" ^
    -db "%DB_NAME%" ^
    -out "%OUTPUT_DIR%\others_vs_morex.out" ^
    -outfmt 6 ^
    -evalue 1e-5


REM ==================================================
REM Step 3: Filter and extract homologs
REM ==================================================
echo.
echo ==================================================
echo Step 3: Filtering and extracting homologs
echo ==================================================

python "%SCRIPT_DIR%\filter_blast_hits.py" ^
    "%OUTPUT_DIR%" ^
    "%INPUT_DIR%\%BARLEY_FASTA%"


REM ==================================================
REM Step 4: Create README
REM ==================================================
echo.
echo ==================================================
echo Step 4: Generating README
echo ==================================================

powershell -ExecutionPolicy Bypass ^
    -File "%SCRIPT_DIR%\generate_readme.ps1" ^
    -OutputDir "%OUTPUT_DIR%"


REM ==================================================
REM Finish
REM ==================================================
echo.
echo ===================================================
echo Pipeline finished successfully!
echo Output directory:
echo %OUTPUT_DIR%
echo ===================================================
pause
exit /b 0


:usage
echo.
echo Usage:
echo   run_thionin_pipeline.bat INPUT_DIR OUTPUT_DIR SCRIPT_DIR
echo.
echo Example:
echo   run_thionin_pipeline.bat input output scripts
echo.
exit /b 1
