@echo off
echo Entpacken von Eigen...
PowerShell -Command "Expand-Archive -Path .\inst\include\Eigen.zip -DestinationPath .\inst\include\Eigen -Force"

echo Entpacken von NLopt...
PowerShell -Command "Expand-Archive -Path .\inst\include\nlopt.zip -DestinationPath .\inst\include\nlopt -Force"

echo Einrichtung abgeschlossen.
pause