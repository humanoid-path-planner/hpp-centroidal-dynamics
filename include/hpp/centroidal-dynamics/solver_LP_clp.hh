/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifdef CLP_FOUND

#ifndef CENTROIDAL_DYNAMICS_LIB_SOLVER_LP_CLP_HH
#define CENTROIDAL_DYNAMICS_LIB_SOLVER_LP_CLP_HH

#include <hpp/centroidal-dynamics/local_config.hh>
#include <hpp/centroidal-dynamics/solver_LP_abstract.hh>
#include <hpp/centroidal-dynamics/util.hh>

#include "ClpSimplex.hpp"

namespace centroidal_dynamics {

class CENTROIDAL_DYNAMICS_DLLAPI Solver_LP_clp : public Solver_LP_abstract {
 private:
  ClpSimplex m_model;

 public:
  Solver_LP_clp();

  /** Solve the linear program
   *  minimize    c' x
   *  subject to  Alb <= A x <= Aub
   *              lb <= x <= ub
   */
  virtual LP_status solve(Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                          Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub,
                          Ref_vectorX sol);

  /** Get the status of the solver. */
  virtual LP_status getStatus();

  /** Get the objective value of the last solved problem. */
  virtual double getObjectiveValue();

  /** Get the value of the dual variables. */
  virtual void getDualSolution(Ref_vectorX res);

  virtual unsigned int getMaximumIterations();

  virtual bool setMaximumIterations(unsigned int maxIter);

  virtual bool setMaximumTime(double seconds);
};

}  // end namespace centroidal_dynamics

#endif  // CENTROIDAL_DYNAMICS_LIB_SOLVER_LP_CLP_HH

#endif  // CLP_FOUND
