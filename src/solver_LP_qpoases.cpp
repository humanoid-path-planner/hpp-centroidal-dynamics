/*
 * Copyright 2015, LAAS-CNRS
 * Author: Andrea Del Prete
 */

#include <robust-equilibrium-lib/solver_LP_qpoases.hh>
#include <robust-equilibrium-lib/logger.hh>

USING_NAMESPACE_QPOASES

namespace robust_equilibrium
{

  Solver_LP_qpoases::Solver_LP_qpoases(): Solver_LP_abstract()
  {
//    m_options.initialStatusBounds = ST_INACTIVE;
//    m_options.setToReliable();
    m_options.setToDefault();
    m_options.printLevel          = PL_NONE; //PL_LOW
    m_options.enableRegularisation = BT_TRUE;
    m_options.enableEqualities = BT_TRUE;
    m_maxIter = 1000;
    m_maxTime = 100.0;
  }

  LP_status Solver_LP_qpoases::solve(Cref_vectorX c, Cref_vectorX lb, Cref_vectorX ub,
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

    int iters = m_maxIter;
    double solutionTime = m_maxTime;
    if(n!=m_solver.getNV() || m!=m_solver.getNC())
    {
      m_solver = SQProblem(n, m, HST_ZERO);
      m_solver.setOptions(m_options);
//      m_solver.printOptions();
      m_init_succeeded = false;
      m_H = MatrixXX::Zero(n,n);
    }

    if(!m_init_succeeded)
    {
      m_status = m_solver.init(NULL, c.data(), A.data(), lb.data(), ub.data(),
                               Alb.data(), Aub.data(), iters, &solutionTime);
      if(m_status==SUCCESSFUL_RETURN)
        m_init_succeeded = true;
    }
    else
    {
      // this doesn't work if I pass NULL instead of m_H.data()
      m_status = m_solver.hotstart(m_H.data(), c.data(), A.data(), lb.data(), ub.data(),
                                  Alb.data(), Aub.data(), iters, &solutionTime);
      if(m_status!=SUCCESSFUL_RETURN)
        m_init_succeeded = false;
    }

    if(m_status==SUCCESSFUL_RETURN)
    {
      m_solver.getPrimalSolution(sol.data());
    }

    return getStatus();
  }

  LP_status Solver_LP_qpoases::getStatus()
  {
    int ss = getSimpleStatus(m_status);
    if(ss==0)
      return LP_STATUS_OPTIMAL;
    if(ss==1)
      return LP_STATUS_MAX_ITER_REACHED;
    if(ss==-2)
       return LP_STATUS_INFEASIBLE;
    if(ss==-3)
      return LP_STATUS_UNBOUNDED;
    return LP_STATUS_ERROR;
  }

  double Solver_LP_qpoases::getObjectiveValue()
  {
    return m_solver.getObjVal();
  }

  void Solver_LP_qpoases::getDualSolution(Ref_vectorX res)
  {
    m_solver.getDualSolution(res.data());
  }

  unsigned int Solver_LP_qpoases::getMaximumIterations()
  {
    return m_maxIter;
  }

  bool Solver_LP_qpoases::setMaximumIterations(unsigned int maxIter)
  {
    if(maxIter==0)
      return false;
    m_maxIter = maxIter;
    return true;
  }

  bool Solver_LP_qpoases::setMaximumTime(double seconds)
  {
    if(seconds<=0.0)
      return false;
    m_maxTime = seconds;
    return true;
  }

} // end namespace robust_equilibrium
