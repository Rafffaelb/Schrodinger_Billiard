@echo off
if exist eigen-3.4.0 (
    echo Eigen library already present.
    exit /b 0
)

echo Downloading Eigen library...
powershell -Command "[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2' -OutFile 'eigen.tar.bz2'"
if %ERRORLEVEL% NEQ 0 (
    echo Failed to download Eigen library.
    exit /b 1
)

tar -xjf eigen.tar.bz2
if %ERRORLEVEL% NEQ 0 (
    echo Failed to extract Eigen library.
    del eigen.tar.bz2
    exit /b 1
)

del eigen.tar.bz2
echo Eigen library downloaded and extracted.