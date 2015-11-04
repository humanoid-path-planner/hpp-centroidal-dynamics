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

You can install cdd lib under Ubuntu 12.04 with the following command:
```
sudo apt-get install libcdd-dev
```
For Eigen and qpOases follow the instructions at the relative webpage (links in the list above).

### Optional
As an alternative to qpOases you can use [CLP](https://projects.coin-or.org/Clp) to solve linear programs.
However, we found qpOases to be faster (especially when solving a series of problems that are similar to each other,
because it can exploit warm start) and more reliable, so we suggest you to stick with qpOases.
In particular, we found that CLP sometimes fails to find the real optimum when using the DLP formulation.
