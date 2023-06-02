# rubicolaMEMS
<p align="center">
A full open-source system that integrates design, simulation and optimization of MEMS sensors.
  <img width=350 length=350 src="https://i.imgur.com/ehGcVeF.png">
</p>

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
sudo add-apt-repository ppa:fenics-packages/fenics
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

### 6. Git clone the repository
It is recommended to make a new directory in your Ubuntu virtual system. Then clone the repository:
```
mkdir rubicola
cd rubicola
git clone https://github.com/ruiesteves-pt/rubicolaMEMS
```

### 7. Installation is done. Run the program (in this case, the simplest accelerometer example) by:
```
python acc_ga.py
```
Make the changes to geometry, figure of merit and FEM studies according to your needs. Further documentation on this in future.

## CITATION
If you use this repository or if this work was useful to you in any way, please consider citing it using the "Cite this repository" tool at the page's the top right or with the following structure:

**MDPI and ACS Style**

Amendoeira Esteves, R.; Wang, C.; Kraft, M. Python-Based Open-Source Electro-Mechanical Co-Optimization System for MEMS Inertial Sensors. Micromachines 2022, 13, 1. https://doi.org/10.3390/mi13010001

**AMA Style**

Amendoeira Esteves R, Wang C, Kraft M. Python-Based Open-Source Electro-Mechanical Co-Optimization System for MEMS Inertial Sensors. Micromachines. 2022; 13(1):1. https://doi.org/10.3390/mi13010001

**Chicago/Turabian Style**

Amendoeira Esteves, Rui, Chen Wang, and Michael Kraft. 2022. "Python-Based Open-Source Electro-Mechanical Co-Optimization System for MEMS Inertial Sensors" Micromachines 13, no. 1: 1. https://doi.org/10.3390/mi13010001 

https://www.mdpi.com/2072-666X/13/1/1

## STATUS
ACTIVE
