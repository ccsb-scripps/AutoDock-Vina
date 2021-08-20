@echo off

set target="Release|x64"
IF "%~1"=="" GOTO noparms
set target="%~1"
:noparms

for /f "usebackq delims=" %%s in (`"%programfiles(x86)%\Microsoft Visual Studio\Installer\vswhere" -latest -property installationPath`) do set devenv="%%s\Common7\IDE\devenv.com"

%devenv% AutoDock-Vina.sln /Build %target%
