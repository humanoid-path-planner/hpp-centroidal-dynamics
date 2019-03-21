/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <hpp/centroidal-dynamics/logger.hh>
#include <hpp/centroidal-dynamics/solver_LP_abstract.hh>
#include <hpp/centroidal-dynamics/solver_LP_qpoases.hh>
#include <iostream>

#ifdef CLP_FOUND
#include <hpp/centroidal-dynamics/solver_LP_clp.hh>
#endif

using namespace std;

namespace centroidal_dynamics {

Solver_LP_abstract *Solver_LP_abstract::getNewSolver(SolverLP solverType) {
  if (solverType == SOLVER_LP_QPOASES) return new Solver_LP_qpoases();

#ifdef CLP_FOUND
  if (solverType == SOLVER_LP_CLP) return new Solver_LP_clp();
#endif

  SEND_ERROR_MSG("Specified solver type not recognized: " + toString(solverType));
  return NULL;
}

bool Solver_LP_abstract::writeLpToFile(const std::string &filename, Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                                       Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub) {
  MatrixXX::Index n = c.size(), m = A.rows();
  assert(lb.size() == n);
  assert(ub.size() == n);
  assert(A.cols() == n);
  assert(Alb.size() == m);
  assert(Aub.size() == m);

  std::ofstream out(filename.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  out.write((const char *)(&n), sizeof(typename MatrixXX::Index));
  out.write((const char *)(&m), sizeof(typename MatrixXX::Index));
  out.write((const char *)c.data(), n * sizeof(typename MatrixXX::Scalar));
  out.write((const char *)lb.data(), n * sizeof(typename MatrixXX::Scalar));
  out.write((const char *)ub.data(), n * sizeof(typename MatrixXX::Scalar));
  out.write((const char *)A.data(), m * n * sizeof(typename MatrixXX::Scalar));
  out.write((const char *)Alb.data(), m * sizeof(typename MatrixXX::Scalar));
  out.write((const char *)Aub.data(), m * sizeof(typename MatrixXX::Scalar));
  out.close();
  return true;
}

bool Solver_LP_abstract::readLpFromFile(const std::string &filename, VectorX &c, VectorX &lb, VectorX &ub, MatrixXX &A,
                                        VectorX &Alb, VectorX &Aub) {
  std::ifstream in(filename.c_str(), std::ios::in | std::ios::binary);
  typename MatrixXX::Index n = 0, m = 0;
  in.read((char *)(&n), sizeof(typename MatrixXX::Index));
  in.read((char *)(&m), sizeof(typename MatrixXX::Index));
  c.resize(n);
  lb.resize(n);
  ub.resize(n);
  A.resize(m, n);
  Alb.resize(m);
  Aub.resize(m);
  in.read((char *)c.data(), n * sizeof(typename MatrixXX::Scalar));
  in.read((char *)lb.data(), n * sizeof(typename MatrixXX::Scalar));
  in.read((char *)ub.data(), n * sizeof(typename MatrixXX::Scalar));
  in.read((char *)A.data(), m * n * sizeof(typename MatrixXX::Scalar));
  in.read((char *)Alb.data(), m * sizeof(typename MatrixXX::Scalar));
  in.read((char *)Aub.data(), m * sizeof(typename MatrixXX::Scalar));
  in.close();
  return true;
}

LP_status Solver_LP_abstract::solve(const std::string &filename, Ref_vectorX sol) {
  VectorX c, lb, ub, Alb, Aub;
  MatrixXX A;
  if (!readLpFromFile(filename, c, lb, ub, A, Alb, Aub)) {
    SEND_ERROR_MSG("Error while reading LP from file " + string(filename));
    return LP_STATUS_ERROR;
  }
  return solve(c, lb, ub, A, Alb, Aub, sol);
}

bool Solver_LP_abstract::setMaximumIterations(unsigned int maxIter) {
  if (maxIter == 0) return false;
  m_maxIter = maxIter;
  return true;
}

bool Solver_LP_abstract::setMaximumTime(double seconds) {
  if (seconds <= 0.0) return false;
  m_maxTime = seconds;
  return true;
}

}  // end namespace centroidal_dynamics
