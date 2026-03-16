@echo off
REM Build acrossCHRphasing with increased stack (64 MB) for Windows
echo Building acrossCHRphasing...
g++ -o acrossCHRphasing.exe acrossCHRphasing.cpp -O2 -fopenmp -lm -Wl,--stack,67108864
if %errorlevel% neq 0 (
    echo Build FAILED
    exit /b 1
)
echo Build OK. Run: acrossCHRphasing.exe pathHap=data/chr pathwindow=data/recombinaisonWindows.txt posOffspring=0 pathresult=./out
