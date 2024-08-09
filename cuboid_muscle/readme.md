# A cuboid muscle model

## Setup
- A dummy cuboid muscle geometry. 
- Solver comprises hyperelasticity solver for the prestretch and both mechanics solver + fastmonodomain solver for the contraction. 
- It uses the CellML model "hodgkin_huxley-razumova".
- No preCICE involved. 

### How to build?
Follow OpenDiHu's documentation for installation, then run 
```
mkorn && sr
```
For a debug build, look into the documentation. 

### How to run?
To run the case go into the build directory and run:

```
./incompressible_mooney_rivlin ../settings_force.py incompressible_mooney_rivlin 10.0
```

> [!WARNING]  
> Currently fails to run in parallel. 
