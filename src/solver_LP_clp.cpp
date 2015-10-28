/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <robust-equilibrium-lib/solver_LP_clp.hh>
#include "CoinBuild.hpp"

namespace robust_equilibrium
{

Solver_LP_clp::Solver_LP_clp(): Solver_LP_abstract()
{
  m_model.setLogLevel(0);
}

LP_status Solver_LP_clp::solve(Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
                               Cref_matrixXX A, Cref_vectorX Alb, Cref_vectorX Aub,
                                Ref_vectorX sol)
{
  int n = (int)c.size();  // number of variables
  int m = (int)A.rows();  // number of constraints
  assert(lb.size()==n);
  assert(ub.size()==n);
  assert(A.cols()==n);
  assert(Alb.size()==m);
  assert(Aub.size()==m);

  m_model.resize(0, n);
  int* rowIndex = new int[n];

  for(int i=0; i<n; i++)
  {
    m_model.setObjectiveCoefficient(i, c(i));
    m_model.setColumnLower(i, lb(i));
    m_model.setColumnUpper(i, ub(i));
    rowIndex[i] = i;
  }

//  m_model.allSlackBasis();

  // This is not the most efficient way to pass the data to the problem
  // but it is the best compromise between efficiency and simplicity.
  // We could be more efficient by using CoinPackedMatrix and ClpPackedMatrix
  // as shown in the example file "addRows.cpp"
  CoinBuild buildObject;
  for (int i=0; i<m; i++)
  {
    buildObject.addRow(n, rowIndex, A.row(i).data(), Alb(i), Aub(i));
  }
  m_model.addRows(buildObject);

  // solve the problem
  m_model.primal();
//  m_model.dual();

  if(m_model.isProvenOptimal())
  {
    const double *solution = m_model.getColSolution();
    for(int i=0; i<n; i++)
      sol(i) = solution[i];
  }

  return getStatus();
}

LP_status Solver_LP_clp::getStatus()
{
  int status = m_model.status();
  if(status<5)
    return (LP_status)status;
  return LP_STATUS_ERROR;
}

double Solver_LP_clp::getObjectiveValue()
{
  return m_model.objectiveValue();
}

unsigned int Solver_LP_clp::getMaximumIterations()
{
  int integerValue;
  m_model.getIntParam(ClpMaxNumIteration, integerValue);
  return integerValue;
}

bool Solver_LP_clp::setMaximumIterations(unsigned int maxIter)
{
  if(maxIter==0)
    return false;
  m_model.setMaximumIterations(maxIter);
  return true;
}

bool Solver_LP_clp::setMaximumTime(double seconds)
{
  if(seconds<=0.0)
    return false;
  m_model.setMaximumSeconds(seconds);
  return true;
}

} // end namespace robust_equilibrium
