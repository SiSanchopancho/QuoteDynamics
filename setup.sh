#!/bin/bash

# Navigate to the directory containing the zip files
cd inst/include

# Unzip Eigen and NLopt
echo "Unzipping Eigen..."
unzip -o Eigen.zip -d Eigen

echo "Unzipping NLopt..."
unzip -o nlopt.zip -d nlopt

# Navigate back to the root directory
cd ../../..

echo "Setup completed."
