@echo off
cd inst\include

echo Unzipping Eigen...
unzip -o Eigen.zip -d Eigen

echo Unzipping NLopt...
unzip -o NLopt.zip -d NLopt

cd ..\..\..

echo Setup completed.
pause
