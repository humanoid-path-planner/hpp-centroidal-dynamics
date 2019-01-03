/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#ifndef CENTROIDAL_DYNAMICS_LIB_SOLVER_LP_ABSTRACT_HH
#define CENTROIDAL_DYNAMICS_LIB_SOLVER_LP_ABSTRACT_HH

#include <Eigen/Dense>
#include <hpp/centroidal-dynamics/config.hh>
#include <hpp/centroidal-dynamics/util.hh>

namespace centroidal_dynamics
{

/**
  * Available LP solvers.
  */
enum CENTROIDAL_DYNAMICS_DLLAPI SolverLP
{
  SOLVER_LP_QPOASES = 0
#ifdef CLP_FOUND
  ,SOLVER_LP_CLP = 1
#endif
};


/**
  * Possible states of an LP solver.
  */
enum CENTROIDAL_DYNAMICS_DLLAPI LP_status
{
  LP_STATUS_UNKNOWN=-1,
  LP_STATUS_OPTIMAL=0,
  LP_STATUS_INFEASIBLE=1,
  LP_STATUS_UNBOUNDED=2,
  LP_STATUS_MAX_ITER_REACHED=3,
  LP_STATUS_ERROR=4
};


/**
 * @brief Abstract interface for a Linear Program (LP) solver.
 */
class CENTROIDAL_DYNAMICS_DLLAPI Solver_LP_abstract
{
protected:
  bool                  m_useWarmStart;   // true if the solver is allowed to warm start
  int                   m_maxIter;        // max number of iterations
  double                m_maxTime;        // max time to solve the LP [s]

public:

  Solver_LP_abstract()
  {
    m_maxIter = 1000;
    m_maxTime = 100.0;
    m_useWarmStart = true;
  }

  /**
   * @brief Create a new LP solver of the specified type.
   * @param solverType Type of LP solver.
   * @return A pointer to the new solver.
   */
  static Solver_LP_abstract* getNewSolver(SolverLP solverType);

  /** Solve the linear program
   *  minimize    c' x
   *  subject to  Alb <= A x <= Aub
   *              lb <= x <= ub
   */
  virtual LP_status solve(Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                          Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub,
                          Ref_vectorX sol) = 0;

  /**
   * @brief Solve the linear program described in the specified file.
   * @param filename Name of the file containing the LP description.
   * @param sol Output solution of the LP.
   * @return A flag describing the final status of the solver.
   */
  virtual LP_status solve(const std::string& filename, Ref_vectorX sol);

  /**
   * @brief Write the specified Linear Program to binary file.
   *  minimize    c' x
   *  subject to  Alb <= A x <= Aub
   *              lb <= x <= ub
   * @param filename
   * @param c
   * @param lb
   * @param ub
   * @param A
   * @param Alb
   * @param Aub
   * @return True if the operation succeeded, false otherwise.
   */
  virtual bool writeLpToFile(const std::string& filename,
                             Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                             Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub);

  /**
   * @brief Read the data describing a Linear Program from the specified binary file.
   * The vectors and matrices are resized inside the method.
   *  minimize    c' x
   *  subject to  Alb <= A x <= Aub
   *              lb <= x <= ub
   * @param filename
   * @param c
   * @param lb
   * @param ub
   * @param A
   * @param Alb
   * @param Aub
   * @return True if the operation succeeded, false otherwise.
   */
  virtual bool readLpFromFile(const std::string& filename,
                              VectorX &c, VectorX &lb, VectorX &ub,
                              MatrixXX &A, VectorX &Alb, VectorX &Aub);

  /** Get the status of the solver. */
  virtual LP_status getStatus() = 0;

  /** Get the objective value of the last solved problem. */
  virtual double getObjectiveValue() = 0;

  /** Get the value of the dual variables. */
  virtual void getDualSolution(Ref_vectorX res) = 0;


  /** Return true if the solver is allowed to warm start, false otherwise. */
  virtual bool getUseWarmStart(){ return m_useWarmStart; }
  /** Specify whether the solver is allowed to use warm-start techniques. */
  virtual void setUseWarmStart(bool useWarmStart){ m_useWarmStart = useWarmStart; }

  /** Get the current maximum number of iterations performed by the solver. */
  virtual unsigned int getMaximumIterations(){ return m_maxIter; }
  /** Set the current maximum number of iterations performed by the solver. */
  virtual bool setMaximumIterations(unsigned int maxIter);


  /** Get the maximum time allowed to solve a problem. */
  virtual double getMaximumTime(){ return m_maxTime; }
  /** Set the maximum time allowed to solve a problem. */
  virtual bool setMaximumTime(double seconds);

};

} // end namespace centroidal_dynamics

#endif //CENTROIDAL_DYNAMICS_LIB_SOLVER_LP_ABSTRACT_HH
