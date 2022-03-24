# How to install preCICE 

The easiest way to install preCICE is by downloading the latest Linux distribution. We used v2.3.0 (libprecice2_2.3.0_focal.deb) - see [preCICE Releases](https://github.com/precice/precice/releases).

First install preCICE system package on Linux / WSL - see [System packages](https://precice.org/installation-packages.html)
Then install the preCICE python bindings globally using the recommended global python installation options:

```
$ pip3 install --user pyprecice
```

Note: make sure you use the same version for the bindings as for the preCICE library.  
Check the installed library using `apt show libprecice2`

If you encounter any python bindings package installation problems, make sure to check [Troubleshooting the pypreCICE installation](https://github.com/precice/python-bindings#troubleshooting--miscellaneous)

It is recommended to then create a local environment in your project folder (/precice_modules for us) and install python bindings and other requirements using:

```
$ cd precice_modules
$ python3 -m venv precice-wsl-env
$ source precice-wsl-env/bin/activate
$ python -m pip install -r requirements.txt
```

There is no need to add the python bindings package pyprecice to the main project requirements as preCICE will not be available in Windows.
See the [aeroelastic_analysis.py](../aeroelastic_analysis.py) example and bash scripts [wsl_script_aero.bash](./wsl_script_aero.bash) and [wsl_script_structures.bash](./wsl_script_structures.bash) for an example of parallelisation of calls to preCICE from Windows.
