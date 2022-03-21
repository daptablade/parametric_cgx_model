First install preCICE system package on Linux / WSL - see (System packages)[https://precice.org/installation-packages.html]
Then install the preCICE python bindings using the recommended global python installation options:

```
pip3 install --user pyprecice
```

Note: make sure you use the same version for the bindings as for the preCICE library.  
Check the installed library using `apt show libprecice2`

If running a python project in Linux, create a local environment in your project folder and install python bindings and other requirements using:

```
python3 -m venv venv
source venv/bin/activate
python -m pip install -r requirements.txt
```

If your python project is windows based, create a local environment, but DO NOT INCLUDE pyprecice in the requirements.
Instead, execute precice using the `wsl -d DISTRO -e python3 ....` command, where DISTRO is the wsl distribution where you installed preCICE. See the `execute_from_windows.py` example for parallelisation of calls to preCICE.

If you encounter any python package installation problems, make sure to check [Troubleshooting the pypreCICE installation](https://github.com/precice/python-bindings#troubleshooting--miscellaneous)
