#!/bin/bash
source precice-wsl-env/bin/activate
wait
echo "Current python version is:" | python --version
echo "Python path is:" | which python

OUTPUT=$(python aeroelastic_two_way.py precice_two_way.xml SolverOne MeshOne)
echo "${OUTPUT}"wsl