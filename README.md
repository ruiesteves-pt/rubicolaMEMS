# rubicolaMEMS
A full 

## Instalation

### 1. Install the Windows Subsystem for Linux (WSL)
Open the Windows Powershell and type in:
```
wsl --install
```
You might have to reboot your computer.
You can verify your installation by opening the command line and type
```
ubuntu
```

### 2. Install Python 3 on your WSL
Open the command line and type ubuntu (this will be required everytime we want to run the program, as it only runs in Ubuntu OS)
Then:
```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install python3
```
Test it by running:
```
python3
```
If you don't ever use Python 2, a handy install would also be:
```
sudo apt install python-is-python3
```
And then test it with:
```
python
```
Then install pip (the package manager):
```
sudo apt install python3-pip
```

### 3. Install FEniCS
To install FEniCS, first add their repository and then install the package:
```
sudo add-apt-repository ppa:fenics-package/fenics
sudo apt-get update
sudo apt-get install fenics
sudo apt-get dist-upgrade
```
You can then test the installation by running with no errors or output:
```
python -c 'import fenics'
```
FEniCS is installed. For more information check: https://github.com/FEniCS/dolfinx#installation

### 4. Install other required libraries
This software requires the use of gmsh, pygmsh, h5py and numpy. Install them using:
```
sudo apt install python3-gmsh
pip install gmsh
pip install pygmsh
pip install numpy
pip install h5py
```

### 5. IF THERE ARE ERRORS
Sometimes some errors arise when using the WSL instead of standalone Ubuntu. The vast majority can be solved by installing the following libraries:
```
sudo apt install mesa-utils libglu1-mesa-dev freeglut3-dev mesa-common-dev
sudo apt-get install libxcursor1
sudo apt-get install libxinerama1
```
