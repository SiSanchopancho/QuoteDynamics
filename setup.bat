@echo off
echo Entpacken von Eigen...
PowerShell -Command "Expand-Archive -Path .\inst\include\Eigen.zip -DestinationPath .\inst\include\Eigen -Force"

echo Entpacken von NLopt...
PowerShell -Command "Expand-Archive -Path .\inst\include\NLopt.zip -DestinationPath .\inst\include\NLopt -Force"

echo Einrichtung abgeschlossen.
pause