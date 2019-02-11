ECHO Storing current directory in VSCMD_START_DIR
SET "VSCMD_START_DIR=%CD%"
CALL "C:\Program Files (x86)\Microsoft Visual Studio\2017\Community\VC\Auxiliary\Build\vcvars64.bat"
cd .build
cmake --build . --config Release
