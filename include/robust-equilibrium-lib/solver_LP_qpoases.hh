/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef ROBUST_EQUILIBRIUM_LIB_SOLVER_QPOASES_HH
#define ROBUST_EQUILIBRIUM_LIB_SOLVER_QPOASES_HH

#include <robust-equilibrium-lib/config.hh>
#include <robust-equilibrium-lib/util.hh>
#include <robust-equilibrium-lib/solver_LP_abstract.hh>
#include <qpOASES.hpp>

namespace robust_equilibrium
{

class ROBUST_EQUILIBRIUM_DLLAPI Solver_LP_qpoases: public Solver_LP_abstract
{
private:
  qpOASES::Options    m_options;  // solver options
  qpOASES::SQProblem  m_solver;   // qpoases solver

  MatrixXX              m_H;              // Hessian matrix
  bool                  m_init_succeeded; // true if solver has been successfully initialized
  qpOASES::returnValue  m_status;         // status code returned by the solver
  int                   m_maxIter;        // max number of iterations
  double                m_maxTime;        // max time to solve the LP [s]

public:

  Solver_LP_qpoases();

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

#endif //ROBUST_EQUILIBRIUM_LIB_SOLVER_QPOASES_HH
