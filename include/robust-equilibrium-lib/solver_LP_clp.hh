/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef ROBUST_EQUILIBRIUM_LIB_SOLVER_LP_CLP_HH
#define ROBUST_EQUILIBRIUM_LIB_SOLVER_LP_CLP_HH

#include <robust-equilibrium-lib/config.hh>
#include <robust-equilibrium-lib/util.hh>
#include <robust-equilibrium-lib/solver_LP_abstract.hh>
#include "ClpSimplex.hpp"

namespace robust_equilibrium
{

class ROBUST_EQUILIBRIUM_DLLAPI Solver_LP_clp: public Solver_LP_abstract
{
private:
  ClpSimplex m_model;

public:

  Solver_LP_clp();

  /** Solve the linear program
   *  minimize    c' x
   *  subject to  Alb <= A x <= Aub
   *              lb <= x <= ub
   */
  LP_status solve(Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                  Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub,
                  Ref_vectorX sol);

  /** Get the status of the solver. */
  LP_status getStatus();

  /** Get the objective value of the last solved problem. */
  double getObjectiveValue();

  void getDualSolution(Ref_vectorX res);

//  void getDualColumnSolution(Ref_vectorX res);

  /** Get the current maximum number of iterations performed
   *  by the solver.
   */
  unsigned int getMaximumIterations();

  /** Set the current maximum number of iterations performed
   *  by the solver.
   */
  bool setMaximumIterations(unsigned int maxIter);

  /** Set the maximum time allowed to solve a problem. */
  bool setMaximumTime(double seconds);

};

} // end namespace robust_equilibrium

#endif //ROBUST_EQUILIBRIUM_LIB_SOLVER_LP_CLP_HH
