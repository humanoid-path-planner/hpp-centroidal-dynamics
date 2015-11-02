/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <robust-equilibrium-lib/solver_LP_abstract.hh>
#include <robust-equilibrium-lib/solver_LP_qpoases.hh>
#include <robust-equilibrium-lib/logger.hh>
#include <iostream>

#ifdef CLP_FOUND
#include <robust-equilibrium-lib/solver_LP_clp.hh>
#endif


using namespace std;

namespace robust_equilibrium
{

Solver_LP_abstract* Solver_LP_abstract::getNewSolver(SolverLP solverType)
{
  if(solverType==SOLVER_LP_QPOASES)
    return new Solver_LP_qpoases();

#ifdef CLP_FOUND
  if(solverType==SOLVER_LP_CLP)
    return new Solver_LP_clp();
#endif

  SEND_ERROR_MSG("Specified solver type not recognized: "+toString(solverType));
  return NULL;
}



} // end namespace robust_equilibrium
