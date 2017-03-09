#~ from centroidal_dynamics import Equilibrium, SolverLP, EquilibriumAlgorithm
from centroidal_dynamics import *

#testing constructors
eq = Equilibrium("test", 54., 4) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES ) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES ) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, False) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, False, 1) 
eq = Equilibrium("test", 54., 4, SolverLP.SOLVER_LP_QPOASES, True, 1, True ) 

#testing all methods
previous = eq.useWarmStart()
eq.setUseWarmStart(False)
assert(previous != eq.useWarmStart())

assert(eq.getName() == "test")


# creating contact points
from numpy import array, asmatrix, matrix

z = array([0.,0.,1.])
P = asmatrix(array([array([x,y,0]) for x in [-0.05,0.05] for y in [-0.1,0.1]]))
N = asmatrix(array([z for _ in range(4)]))

eq.setNewContacts(asmatrix(P),asmatrix(N),0.3,EquilibriumAlgorithm.EQUILIBRIUM_ALGORITHM_LP)

c= asmatrix(array([0.,0.,1.]))
ddc= asmatrix(array([0.,0.,0.]))

status, robustness = eq.computeEquilibriumRobustness(c,ddc)
assert (status == LP_STATUS_OPTIMAL), "LP should not fail"
assert (robustness > 0), "first test should be in equilibrirum"
	
ddc= asmatrix(array([1000.,0.,0.]))
status, robustness = eq.computeEquilibriumRobustness(c,ddc)
assert (status == LP_STATUS_OPTIMAL), "LP should not fail"
assert (robustness < 0), "first test should NOT be in equilibrirum"
