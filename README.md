## Introduction

This is my implementation of FEM cloth simulation assignment in [CSC417/CSC2549-Physics-based Animation](https://github.com/dilevin/CSC417-physics-based-animation).
![Cloth simulation!](images/cloth.gif)

## Build & Execution
```
git submodule update --init --recursive
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make 
./a4-cloth-simulation
```