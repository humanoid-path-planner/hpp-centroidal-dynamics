# RobustEquilibriumLib
Utility classes to check the (robust) equilibrium of a system in contact with the environment.
The main class that collects all the equilibrium-related algorithms is ```StaticEquilibrium```.
All the algorithms take as input:
* A list of contact points
* A list of contact normals
* The contact friction coefficient
* The number of generators used for the linear approximations of the friction cones
* The mass of the system

Once these input parameters have been specified, the user has access to four algorithms implemented in the following four methods of the class ```StaticEquilibrium```:
* ```computeEquilibriumRobustness```: compute the robustness of the equilibrium of a given a CoM (center of mass) position (negative values mean the system can not be in equilibrium).
* ```checkRobustEquilibrium```: checks whether the system can be in equilibrium in a specified CoM position and with a specified robustness level.
* ```findExtremumOverLine```: Find the extremum com position that is in robust equilibrium along the specified line.
* ```findExtremumInDirection```: Find the extremum com position that is in robust equilibrium in the specified direction.

All these problems boil down to solving Linear Programs.
Different formulations are implemented and tested in ```test_static_equilibrium```.
More details can be found in the code documentation.
In the end, we found that most of the times the dual LP formulation (DLP) is the fastest.

The test ```test_LP_solvers``` tries to solve some LP problems using qpOases and checks that the results are correct.

## Dependencies
* [Eigen (version >= 3.2.2)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
* [cdd lib](https://www.inf.ethz.ch/personal/fukudak/cdd_home/)
* [qpOases (version >= 3.0beta)](https://projects.coin-or.org/qpOASES)

## Installation Steps for Ubuntu 12.04
You can install cdd lib with the following command:
```
sudo apt-get install libcdd-dev
```
You can install Eigen3 with the following command:
```
 sudo apt-get install libeigen3-dev
```
For qpOases you have to install the pkg-config version you can find here: https://github.com/humanoid-path-planner/qpoases
```
git clone --recursive https://github.com/humanoid-path-planner/qpoases
mkdir qpoases/build
cd qpoases/build
cmake -DCMAKE_INSTALL_PREFIX=${DEVEL_DIR}/install ..
make install
```
Then you can clone this repository using ssh:
```
git clone --recursive git@github.com:andreadelprete/robust-equilibrium-lib.git $ROBUST_EQUI_LIB_DIR
```
or using http:
```
git clone --recursive https://github.com/andreadelprete/robust-equilibrium-lib.git $ROBUST_EQUI_LIB_DIR
```
And you can build this library using CMake:
```
mkdir $ROBUST_EQUI_LIB_DIR/build
cd $ROBUST_EQUI_LIB_DIR/build
cmake -DCMAKE_INSTALL_PREFIX=${DEVEL_DIR}/install ..
make install
```
Currently, CMake may have problems finding CDD.
If this is the case you can specify its path manually, for instance:
```
cmake -DCDD_LIBRARY=/usr/lib/libcdd.so -DCMAKE_INSTALL_PREFIX=${DEVEL_DIR}/install ..
```

### Optional
As an alternative to qpOases you can use [CLP](https://projects.coin-or.org/Clp) to solve linear programs.
However, we found qpOases to be faster (especially when solving a series of problems that are similar to each other,
because it can exploit warm start) and more reliable, so we suggest you to stick with qpOases.
In particular, we found that CLP sometimes fails to find the real optimum when using the DLP formulation.
