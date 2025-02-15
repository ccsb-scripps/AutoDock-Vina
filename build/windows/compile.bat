@echo off
setlocal enabledelayedexpansion

:: Default settings
set config=Release
set platform=x64

:: Default location of MSBuild.exe
set "msbuild_path=C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\MSBuild\Current\Bin\amd64\MSBuild.exe"

:: Allow user to override with command-line arguments
if NOT "%~1"=="" set config=%~1
if NOT "%~2"=="" set platform=%~2

:: Check if Boost library is provided
if NOT "%~3"=="" (
    set "boost_lib=%~3"
) else (
    :: Default Boost library path based on choice of platform
    if /i "%platform%"=="x64" (
        set boost_lib=C:\local\boost_1_83_0\lib64-msvc-14.3
    ) else (
        set boost_lib=C:\local\boost_1_83_0\lib32-msvc-14.3
    )
)

:: Validate Boost library path
if NOT exist "%boost_lib%" (
    echo ERROR: Boost library files not found at "%boost_lib%" for %platform%!
    exit /b 1
)
echo Using Boost library: "%boost_lib%"

:: Check if MSBuild path is provided
if NOT "%~4"=="" (
    set "msbuild_path=%~4"
)

:: Validate MSBuild path
if NOT exist "%msbuild_path%" (
    echo WARNING: MSBuild not found at "%msbuild_path%". Searching...

    set msbuild_path=

    :: Define possible search locations
    set search_dirs="%ProgramFiles(x86)%" "%ProgramFiles%" "%LocalAppData%"
    set expect_paths=MSBuild\Current\Bin\amd64\MSBuild.exe

    if /i "%platform%" NEQ "x64" set expect_paths=MSBuild\Current\Bin\MSBuild.exe

    :: Loop through search directories to find MSBuild
    for %%D in (%search_dirs%) do (
        for /f "delims=" %%A in ('cmd /c dir "%%D\Microsoft Visual Studio\" /s /b 2^>nul') do (
            if exist "%%A\%expect_paths%" (
                set "msbuild_path=%%A\%expect_paths%"
                goto found
            )
        )
    )

    :: If MSBuild.exe is not found
    if not defined msbuild_path (
        echo ERROR: MSBuild.exe not found for %platform%!
        exit /b 1
    )
)

:: If MSBuild.exe is found
:found
echo Using MSBuild: "%msbuild_path%"

:: Retrieve Git Version
set GIT_VERSION=0.0.0-unknown
for /f %%i in ('git describe --abbrev=7 --dirty --always --tags 2^>nul') do set GIT_VERSION=%%i

:: Print Git Version
echo Program Version from Git: "%GIT_VERSION%"

:: Run MSBuild
call "%msbuild_path%" AutoDock-Vina.sln /p:Configuration=%config% /p:Platform=%platform% /p:BoostLibraryPath="%boost_lib%" /p:GIT_VERSION="%GIT_VERSION%" /m
endlocal
